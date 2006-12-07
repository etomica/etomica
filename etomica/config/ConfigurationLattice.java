package etomica.config;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorArrayListCompound;
import etomica.integrator.IntegratorHard;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorSequential;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.SpaceLattice;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Constructs configuration that has the molecules placed on the sites of a
 * lattice. Molecules are placed sequentially on the sites of the lattice in an
 * order specified by an index iterator set at construction.
 * <p>
 * Iteration over the indexes yields integer arrays, and with each iteration the
 * array is passed to the site method of the lattice which returns the position
 * for placement of the next molecule. Each array index is iterated to a maximum
 * value determined by the number of molecules to be placed, the dimensions of
 * the phase in which they are placed, and the lattice constants of the lattice.
 * <p>
 * An instance of this class may be configured to place atoms such that they
 * uniformly fill the volume of the phase. It will attempt this by scaling the
 * lattice constants of the configuration in an appropriate way. Success in
 * getting a good spatial distribution may vary.
 * <p>
 * An instance can also be configured to remember the indices used to get
 * lattice position for each molecule that is placed. This can be useful if it
 * is desired to associate each molecule with a lattice site.
 */
public class ConfigurationLattice extends Configuration {

    /**
     * Constructs class using instance of IndexIteratorSequential as the default
     * index iterator.
     */
    public ConfigurationLattice(SpaceLattice lattice) {
        this(lattice, new IndexIteratorSequential(lattice.D()));
    }

    /**
     * Construct class that will place atoms on sites of the given lattice,
     * proceeding in the order resulting from iteration through the given index
     * iterator.
     */
    public ConfigurationLattice(SpaceLattice lattice,
            IndexIteratorSizable indexIterator) {
        super(lattice.getSpace());
        if(indexIterator.getD() != lattice.D()) {
            throw new IllegalArgumentException("Dimension of index iterator and lattice are incompatible");
        }
        this.lattice = lattice;
        this.indexIterator = indexIterator;
        atomActionTranslateTo = new AtomActionTranslateTo(lattice.getSpace());
        atomIterator = new AtomIteratorArrayListCompound();
        work = space.makeVector();
    }

    /**
     * Specifies whether indices used to place each atom should be remembered.
     * Default is false.
     */
    public void setRememberingIndices(boolean b) {
        indices = b ? new int[0][] : null;
    }

    /**
     * Indicates of instance is set to remember indices used to place each atom.
     */
    public boolean isRememberingIndices() {
        return indices != null;
    }

    /**
     * Returns an array of integer arrays, each corresponding to the indexes
     * used to place a molecule in the last call to initializePositions. The
     * index for each atom is obtained from the returned array using the atom's
     * global index; that is, <code>getIndices()[atom.getGlobalIndex()]</code>
     * gives the index array for <code>atom</code>.  A new instance of this
     * array is made with each call to initializeCoordinates.
     * 
     * @throws IllegalStateException
     *             if isRememberingIndices if false.
     */
    public int[][] getIndices() {
        if (indices == null) {
            throw new IllegalStateException(
                    "ConfigurationLattice is not set to remember the indices.");
        }
        return indices;
    }

    /**
     * Places the molecules in the given phase on the positions of the
     * lattice.  
     */
    public void initializeCoordinates(Phase p) {
        if (indices != null) {
            indices = new int[p.getSpeciesMaster().getMaxGlobalIndex()+1][];
        }
        super.initializeCoordinates(p);
    }

    public void initializePositions(AtomArrayList[] lists) {
        if (lists.length == 0)
            return;
        atomIterator.setLists(lists);
        int sumOfMolecules = atomIterator.size();
        if (sumOfMolecules == 0) {
            return;
        }
        if (sumOfMolecules == 1) {
            atomIterator.reset();
            work.E(0.0);
            atomActionTranslateTo.setDestination(work);
            atomActionTranslateTo.actionPerformed(atomIterator.nextAtom());
            return;
        }
        int basisSize = 1;
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int) Math.ceil((double) sumOfMolecules
                / (double) basisSize);

        // determine scaled shape of simulation volume
        double[] shape = (double[]) dimensions.clone();
        double[] latticeConstant = lattice.getLatticeConstants();
        for (int i = 0; i < shape.length; i++) {
            shape[i] /= latticeConstant[i];
        }

        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        if (indexIterator.getD() > latticeDimensions.length) {
            int[] iteratorDimensions = new int[latticeDimensions.length+1];
            System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0,
                    latticeDimensions.length);
            iteratorDimensions[latticeDimensions.length] = basisSize;
            indexIterator.setSize(iteratorDimensions);
        }
        else {
            indexIterator.setSize(latticeDimensions);
        }

        // determine lattice constant
        Vector latticeScaling = space.makeVector();
        if (rescalingToFitVolume) {
            for (int i = 0; i < latticeDimensions.length; i++) {
                // in favorable situations, this should be approximately equal
                // to 1.0
                latticeScaling.setX(i, dimensions[i]
                        / (latticeDimensions[i] * latticeConstant[i]));
            }
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        double[] offset = (double[]) dimensions.clone();
        double[] vectorOfMax = new double[lattice.getSpace().D()];
        double[] vectorOfMin = new double[lattice.getSpace().D()];
        for (int i = 0; i < lattice.getSpace().D(); i++) {
            vectorOfMax[i] = Double.NEGATIVE_INFINITY;
            vectorOfMin[i] = Double.POSITIVE_INFINITY;
        }

        // XXX this looks scary and asking for trouble
        // it's probably not needed/wanted for periodic boundaries, but
        // gets the atoms off the boundaries for non-periodic boundaries
        indexIterator.reset();
        while (indexIterator.hasNext()) {
            Vector site = (Vector) lattice.site(indexIterator.next());
            site.TE(latticeScaling);
            for (int i = 0; i < site.D(); i++) {
                vectorOfMax[i] = Math.max(vectorOfMax[i], site.x(i));
                vectorOfMin[i] = Math.min(vectorOfMin[i], site.x(i));
            }
        }
        for (int i = 0; i < lattice.getSpace().D(); i++) {
            offset[i] = -0.5 * (vectorOfMax[i] - vectorOfMin[i])
                    - vectorOfMin[i];
        }
        Vector offsetVector = Space.makeVector(offset);

        myLat = new MyLattice(lattice, latticeScaling, offsetVector);

        // Place molecules
        atomIterator.reset();
        indexIterator.reset();
        while (atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            if (!a.getNode().isLeaf()) {
                // initialize coordinates of child atoms
                Conformation config = a.getType().creator().getConformation();
                config.initializePositions(((AtomTreeNodeGroup) a.getNode()).getChildList());
            }

            int[] ii = indexIterator.next();
            Vector site = (Vector) myLat.site(ii);
            if (indices != null) {
                indices[a.getGlobalIndex()] = (int[]) ii.clone();
            }
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
    }

    private int[] calculateLatticeDimensions(int nCells, double[] shape) {
        int dimLeft = shape.length;
        int nCellsLeft = nCells;
        int[] latticeDimensions = new int[shape.length];
        while (dimLeft > 0) {
            double smin = Double.POSITIVE_INFINITY;
            int dmin = 0;
            double product = 1.0;
            for (int idim = 0; idim < shape.length; idim++) {
                if (latticeDimensions[idim] > 0)
                    continue;
                if (shape[idim] < smin) {
                    smin = shape[idim];
                    dmin = idim;
                }
                product *= shape[idim];
            }
            // round off except for last dimension (then round up)
            if (dimLeft > 1) {
                latticeDimensions[dmin] = (int) Math.round(shape[dmin]
                        * Math.pow((nCellsLeft / product), 1.0 / dimLeft));
            } else {
                latticeDimensions[dmin] = (int) Math.ceil(shape[dmin]
                        * nCellsLeft / product);
            }
            nCellsLeft = (nCellsLeft + latticeDimensions[dmin] - 1)
                    / latticeDimensions[dmin];
            dimLeft--;
        }
        return latticeDimensions;
    }

    /**
     * Returns a lattice with positions the same as those used in the 
     * most recent use of initializeCoordinates.  Includes any scaling
     * or translation applied to fill the space, and thus will not necessarily
     * be the same positions as specified by the lattice given at construction.
     */
    public SpaceLattice getLatticeMemento() {
        return myLat;
    }

    private final SpaceLattice lattice;
    private final IndexIteratorSizable indexIterator;
    private final Vector work;
    private boolean rescalingToFitVolume = true;
    private final AtomActionTranslateTo atomActionTranslateTo;
    private final AtomIteratorArrayListCompound atomIterator;
    private MyLattice myLat;
    private int[][] indices = null;
    private static final long serialVersionUID = 1L;

    // /**
    // * Sets the size of the lattice (number of atoms in each direction) so
    // that
    // * it has a number of sites equal or greater than the given value n. Sets
    // * lattice to be square (same number in all directions), and finds
    // smallest
    // * size that gives number of sites equal or exceeding the given value.
    // */
    // public void setSize(int n) {
    // if (n < 0 || n >= Integer.MAX_VALUE)
    // throw new IllegalArgumentException(
    // "Inappropriate size specified for lattice: " + n);
    // int i = (int)Math.pow(n,1/D());
    // while(i <= n) {
    // if ((int)Math.pow(i,D()) >= n) break;
    // i++;
    // }
    // setSize(primitive.space.makeArrayD(i));
    // }
    //
    // /**
    // * Performs same actions as setSize(int), then size of primitive vectors
    // are
    // * adjusted such that lattice will fit in the given dimensions. Assumes
    // all
    // * dimension values are equal (future development will address case of
    // non-
    // * square dimensions).
    // */
    // public void setSize(int n, double[] dimensions) {
    // if (n < 0 || n >= Integer.MAX_VALUE)
    // throw new IllegalArgumentException(
    // "Inappropriate size specified for lattice: " + n);
    // int i = (int)Math.pow(n,1/D());
    // while(i <= n) {
    // if ((int)Math.pow(i,D()) >= n) break;
    // i++;
    // }
    // //size primitive vectors to given dimensions
    // double[] newLength = new double[D()];
    // for (int j = 0; j < D(); j++)
    // newLength[j] = dimensions[j] / (double) i;
    // primitive.setSize(newLength);
    // setSize(primitive.space.makeArrayD(i));
    // }

    public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
        sim.getDefaults().atomSize = 5.0;
        Phase phase = new Phase(sim);
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
        int k = 4;
        phase.getAgent(species).setNMolecules(4 * k * k * k);
        IntegratorHard integrator = new IntegratorHard(sim);
        integrator.setPhase(phase);
//        ColorSchemeByType colorScheme = new ColorSchemeByType();
        // CubicLattice lattice = new LatticeCubicBcc();
        BravaisLatticeCrystal lattice = new LatticeCubicFcc();
        // CubicLattice lattice = new LatticeCubicSimple();
        ConfigurationLattice configuration = new ConfigurationLattice(lattice);
        // phase.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
        configuration.initializeCoordinates(phase);
        // etomica.graphics.DisplayPhase display = new
        // etomica.graphics.DisplayPhase(phase);

        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(
                sim);
//        ((ColorSchemeByType) ((DisplayPhase) simGraphic.displayList()
//                .getFirst()).getColorScheme()).setColor(species
//                .getMoleculeType(), java.awt.Color.red);
        simGraphic.makeAndDisplayFrame();
    }

    /**
     * Returns the resizeLatticeToFitVolume flag.
     * 
     * @return boolean
     */
    public boolean isRescalingToFitVolume() {
        return rescalingToFitVolume;
    }

    /**
     * Sets the resizeLatticeToFitVolume flag, which if true indicates that the
     * primitive vectors should be resized to fit the dimensions of the phase.
     * Default is true.
     * 
     * @param resizeLatticeToFitVolume
     *            The resizeLatticeToFitVolume to set
     */
    public void setRescalingToFitVolume(boolean resizeLatticeToFitVolume) {
        this.rescalingToFitVolume = resizeLatticeToFitVolume;
    }

    /**
     * Used to store the state of a lattice.
     * 
     * @author nancycribbin, Andrew Schultz, Dr. Kofke
     * 
     */
    public static class MyLattice implements SpaceLattice {

        MyLattice(SpaceLattice l, Vector latticeScaling, Vector offset) {
            lattice = l;
            this.latticeScaling = (Vector) latticeScaling.clone();
            this.offset = (Vector) offset.clone();

        }

        public Space getSpace() {
            return lattice.getSpace();
        }

        public int D() {
            return lattice.D();
        }

        public Object site(int[] index) {
            Vector site = (Vector) lattice.site(index);
            site.TE(latticeScaling);
            site.PE(offset);

            return site;
        }

        public double[] getLatticeConstants() {
            double[] lat = lattice.getLatticeConstants();
            for (int i = 0; i < lat.length; i++) {
                lat[i] *= latticeScaling.x(i);
            }
            return lat;
        }

        SpaceLattice lattice;
        public Vector latticeScaling;
        Vector offset;

    }

}
