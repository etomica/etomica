package etomica.config;

import etomica.action.AtomActionTranslateTo;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IVector;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.integrator.IntegratorHard;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.SpaceLattice;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
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
 * the box in which they are placed, and the lattice constants of the lattice.
 * If there are not enough molecules to fill the lattice, they will be
 * distributed as evenly as possible.
 * <p>
 * An instance of this class may be configured to place atoms such that they
 * uniformly fill the volume of the box. It will attempt this by scaling the
 * lattice constants of the configuration in an appropriate way. Success in
 * getting a good spatial distribution may vary.
 * <p>
 * An instance can also be configured to remember the indices used to get
 * lattice position for each molecule that is placed. This can be useful if it
 * is desired to associate each molecule with a lattice site.
 */
public class ConfigurationLattice implements Configuration, java.io.Serializable {

    /**
     * Constructs class using instance of IndexIteratorRectangular as the default
     * index iterator.
     */
    public ConfigurationLattice(SpaceLattice lattice, Space space) {
        this(lattice, new IndexIteratorRectangular(lattice.D()), space);
    }

    /**
     * Construct class that will place atoms on sites of the given lattice,
     * proceeding in the order resulting from iteration through the given index
     * iterator.
     */
    public ConfigurationLattice(SpaceLattice lattice,
            IndexIteratorSizable indexIterator, Space space) {
        if(indexIterator.getD() != lattice.D()) {
            throw new IllegalArgumentException("Dimension of index iterator and lattice are incompatible");
        }
        this.lattice = lattice;
        this.indexIterator = indexIterator;
        this.space = space;
        atomActionTranslateTo = new AtomActionTranslateTo(lattice.getSpace());
        setBoundaryPadding(0);
    }

    public void setBoundaryPadding(double newBoundaryPadding) {
        boundaryPadding = newBoundaryPadding;
    }
    
    public double getBoundaryPadding() {
        return boundaryPadding;
    }

    /**
     * Places the molecules in the given box on the positions of the
     * lattice.  
     */
    public void initializeCoordinates(IBox box) {
        IAtomSet moleculeList = box.getMoleculeList();
        int sumOfMolecules = moleculeList.getAtomCount();
        if (sumOfMolecules == 0) {
            return;
        }
        int basisSize = 1;
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int) Math.ceil((double) sumOfMolecules
                / (double) basisSize);

        // determine scaled shape of simulation volume
        IVector shape = space.makeVector();
        shape.E(box.getBoundary().getDimensions());
        shape.PE(-boundaryPadding);
        IVector latticeConstantV = Space.makeVector(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        int nSites = basisSize;
        for (int i=0; i<latticeDimensions.length; i++) {
            nSites *= latticeDimensions[i];
        }
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
        IVector latticeScaling = space.makeVector();
        if (rescalingToFitVolume) {
            // in favorable situations, this should be approximately equal
            // to 1.0
            latticeScaling.E(box.getBoundary().getDimensions());
            latticeScaling.PE(-boundaryPadding);
            latticeScaling.DE(latticeConstantV);
            latticeScaling.DE(Space.makeVector(latticeDimensions));
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        IVector offset = space.makeVector();
        offset.E(box.getBoundary().getDimensions());
        IVector vectorOfMax = space.makeVector();
        IVector vectorOfMin = space.makeVector();
        IVector site = space.makeVector();
        vectorOfMax.E(Double.NEGATIVE_INFINITY);
        vectorOfMin.E(Double.POSITIVE_INFINITY);

        // XXX this can do strange things. it's probably not needed for 
        // periodic boundaries, but gets the atoms off the boundaries for 
        // non-periodic boundaries
        indexIterator.reset();

        while (indexIterator.hasNext()) {
            site.E((IVector) lattice.site(indexIterator.next()));
            site.TE(latticeScaling);
            for (int i=0; i<site.getD(); i++) {
                vectorOfMax.setX(i, Math.max(site.x(i),vectorOfMax.x(i)));
                vectorOfMin.setX(i, Math.min(site.x(i),vectorOfMin.x(i)));
            }
        }
        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
        offset.TE(-0.5);
        offset.ME(vectorOfMin);

        myLat = new MyLattice(lattice, latticeScaling, offset);

        // Place molecules
        indexIterator.reset();
        double voidFrac = (nSites - sumOfMolecules)/((double)nSites);
        double voidSum = 0;
        int siteCount = 0;
        for (int i=0; i<sumOfMolecules; i++) {
            IMolecule a = (IMolecule)moleculeList.getAtom(i);
            int[] ii = indexIterator.next();
            siteCount++;
            // add voidFrac for each /site/ (not molecule)
            voidSum += voidFrac;
            while (voidSum > 1.0) {
                // we've gone through enough sites that we should insert a void
                // now.  Subtract one, but still add voidFrac since we're still
                // advancing one site.
                voidSum += voidFrac - 1;
                ii = indexIterator.next();
                siteCount++;
            }
            // initialize coordinates of child atoms
        	atomActionTranslateTo.setAtomPositionDefinition(a.getType().getPositionDefinition());
            Conformation config = ((AtomTypeMolecule)a.getType()).getConformation();
            config.initializePositions(a.getChildList());

            atomActionTranslateTo.setDestination((IVector)myLat.site(ii));
            atomActionTranslateTo.actionPerformed(a);
        }
        if (nSites - siteCount > 1) {
            // nSites - siteCount = 0 is ideal.
            // indexIterator.next() would throw if nSites < siteCount
            // nSites - siteCount = 1 will be typical for cases where the void distribution can't be perfect
            // so we just need to check for nSites - siteCount > 1
            throw new RuntimeException("Failed to properly iterate through the lattice sites "+nSites+" "+siteCount);
        }
    }

    protected int[] calculateLatticeDimensions(int nCells, IVector shape) {
        int dimLeft = shape.getD();
        int nCellsLeft = nCells;
        int[] latticeDimensions = new int[shape.getD()];
        while (dimLeft > 0) {
            double smin = Double.POSITIVE_INFINITY;
            int dmin = 0;
            double product = 1.0;
            for (int idim = 0; idim < shape.getD(); idim++) {
                if (latticeDimensions[idim] > 0)
                    continue;
                if (shape.x(idim) < smin) {
                    smin = shape.x(idim);
                    dmin = idim;
                }
                product *= shape.x(idim);
            }
            // round off except for last dimension (then round up)
            if (dimLeft > 1) {
                latticeDimensions[dmin] = (int) Math.round(shape.x(dmin)
                        * Math.pow((nCellsLeft / product), 1.0 / dimLeft));
            } else {
                latticeDimensions[dmin] = (int) Math.ceil(shape.x(dmin)
                        * nCellsLeft / product);
            }
            if (latticeDimensions[dmin] == 0){
            	latticeDimensions[dmin] = 1;
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

    protected final SpaceLattice lattice;
    protected final IndexIteratorSizable indexIterator;
    protected boolean rescalingToFitVolume = true;
    protected final AtomActionTranslateTo atomActionTranslateTo;
    protected MyLattice myLat;
    protected double boundaryPadding;
    protected final Space space;
    private static final long serialVersionUID = 3L;

    public static void main(String[] args) {
    	Space sp = Space3D.getInstance();
        Simulation sim = new Simulation(sp);
        IPotentialMaster potentialMaster = new PotentialMaster(sim.getSpace());
        IBox box = new Box(sim, sp);
        sim.addBox(box);
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species);
        ((AtomTypeSphere)species.getMoleculeType().getChildTypes()[0]).setDiameter(5.0);
        int k = 4;
        box.setNMolecules(species, 4 * k * k * k);
        IntegratorHard integrator = new IntegratorHard(sim, potentialMaster, sp);
        integrator.setBox(box);
//        ColorSchemeByType colorScheme = new ColorSchemeByType();
        // CubicLattice lattice = new LatticeCubicBcc();
        BravaisLatticeCrystal lattice = new LatticeCubicFcc();
        // CubicLattice lattice = new LatticeCubicSimple();
        ConfigurationLattice configuration = new ConfigurationLattice(lattice, sp);
        // box.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
        configuration.initializeCoordinates(box);
        // etomica.graphics.DisplayBox display = new
        // etomica.graphics.DisplayBox(box);

        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(
                sim, sp);
//        ((ColorSchemeByType) ((DisplayBox) simGraphic.displayList()
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
     * primitive vectors should be resized to fit the dimensions of the box.
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

        public MyLattice(SpaceLattice l, IVector latticeScaling, IVector offset) {
            lattice = l;
            this.latticeScaling = latticeScaling;
            this.offset = offset;
            this.site = l.getSpace().makeVector();
        }

        public Space getSpace() {
            return lattice.getSpace();
        }

        public int D() {
            return lattice.D();
        }

        /**
         * Returns the same instance of IVector with each call.
         */
        public Object site(int[] index) {
            site.E((IVector) lattice.site(index));
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

        final SpaceLattice lattice;
        final public IVector latticeScaling;
        final IVector offset;
        final IVector site;

    }

}
