package etomica.paracetamol;

import etomica.action.AtomActionTranslateTo;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.config.Configuration;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.SpaceLattice;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
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
public class ConfigurationMonoclinicLattice extends Configuration implements AgentSource {

	private final static String APP_NAME = "Configuration Monoclinic Lattice";

    /**
     * Constructs class using instance of IndexIteratorRectangular as the default
     * index iterator.
     */
    public ConfigurationMonoclinicLattice(SpaceLattice lattice) {
        this(lattice, new IndexIteratorRectangular(lattice.D()));
    }

    /**
     * Construct class that will place atoms on sites of the given lattice,
     * proceeding in the order resulting from iteration through the given index
     * iterator.
     */
    public ConfigurationMonoclinicLattice(SpaceLattice lattice,
            IndexIteratorSizable indexIterator) {
        if(indexIterator.getD() != lattice.D()) {
            throw new IllegalArgumentException("Dimension of index iterator and lattice are incompatible");
        }
        this.lattice = lattice;
        this.indexIterator = indexIterator;
        atomActionTranslateTo = new AtomActionTranslateTo(lattice.getSpace());
        atomGroupAction = new AtomGroupAction(new AtomActionTransformed(lattice.getSpace()));
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
     *             if isRememberingIndices if false, or if initializeCoordinates has not previously been called.
     */
    public int[][] getIndices() {
        if (indices == null) {
            throw new IllegalStateException(
                    "ConfigurationLattice is not set to remember the indices.");
        }
        if (indices[0].length == 0) {
            throw new IllegalStateException("initializeCoordinates has not been called for this phase");
        }
        return indices;
    }

    public Class getAgentClass() {
		return IVector.class;
	}

	public Object makeAgent(IAtom a) {
		return null;
	}

	public void releaseAgent(Object agent, IAtom atom) {
		
	}

	public AtomAgentManager getCellManager(){
		
		return cellManager;
	}
    
    /**
     * Places the molecules in the given phase on the positions of the
     * lattice.  
     */
    public void initializeCoordinates(Phase phase) {
        if (indices != null) {
            indices = new int[phase.getSpeciesMaster().getMaxGlobalIndex()+1][];
        }
        AtomIteratorAllMolecules atomIterator = new AtomIteratorAllMolecules(phase);
        int sumOfMolecules = atomIterator.size();
        if (sumOfMolecules == 0) {
            return;
        }
        int basisSize = 4;
        
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int) Math.ceil((double) sumOfMolecules
                / (double) basisSize);

        // determine scaled shape of simulation volume
        IVector shape = phase.getSpace().makeVector();
        shape.E(phase.getBoundary().getDimensions());
        IVector latticeConstantV = Space.makeVector(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

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
        IVector latticeScaling = phase.getSpace().makeVector();
        if (rescalingToFitVolume) {
            // in favorable situations, this should be approximately equal
            // to 1.0
            latticeScaling.E(phase.getBoundary().getDimensions());
            latticeScaling.DE(latticeConstantV);
            latticeScaling.DE(Space.makeVector(latticeDimensions));
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        IVector offset = phase.getSpace().makeVector();
        offset.E(phase.getBoundary().getDimensions());
        IVector vectorOfMax = phase.getSpace().makeVector();
        IVector vectorOfMin = phase.getSpace().makeVector();
        IVector site = phase.getSpace().makeVector();
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
        atomIterator.reset();
        indexIterator.reset();

    	ConformationParacetamolMonoclinic regConfig = new ConformationParacetamolMonoclinic(lattice.getSpace());
    	cellManager = new AtomAgentManager(this, phase);
    	IVector cellPosition = null;
    	Tensor t = lattice.getSpace().makeTensor();
    	
    	for (IAtomGroup atom = (IAtomGroup)atomIterator.nextAtom(); atom != null;
        	atom = (IAtomGroup)atomIterator.nextAtom()) {
    		    
            int[] ii = indexIterator.next();
            
            regConfig.initializePositions(atom.getChildList());
            
            switch (ii[3]){
            case 0:
            	t.setComponent(0, 0, 1);
            	t.setComponent(1, 1, 1);
            	t.setComponent(2, 2, 1);
            	break;
            case 1:
            	t.setComponent(0, 0,-1);
            	t.setComponent(1, 1, 1);
            	t.setComponent(2, 2,-1);
            	break;
            case 2:
            	t.setComponent(0, 0,-1);
            	t.setComponent(1, 1,-1);
            	t.setComponent(2, 2,-1);
            	break;
            case 3:
            	t.setComponent(0, 0, 1);
            	t.setComponent(1, 1,-1);
            	t.setComponent(2, 2, 1);
            	break;
            }
      
      	  ((AtomActionTransformed)atomGroupAction.getAction()).setTransformationTensor(t);
          atomGroupAction.actionPerformed(atom);
 
            if (indices != null) {
                indices[atom.getGlobalIndex()] = (int[]) ii.clone();
            }
            atomActionTranslateTo.setDestination((IVector)myLat.site(ii));
            atomActionTranslateTo.actionPerformed(atom);
            
            if (ii[3] == 0){
            	cellPosition = phase.getSpace().makeVector();
            	cellPosition.E((IVector)myLat.site(ii));
            	
            	
            //remember the coordinate of the cell
            //Loop 8 times over the basis and we can make the cell assignment here!!
            }
            cellManager.setAgent(atom, cellPosition);
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
                if (latticeDimensions[dmin] == 0){
                	latticeDimensions[dmin] = 1;
                }
            } else {
                latticeDimensions[dmin] = (int) Math.ceil(shape.x(dmin)
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

    protected final SpaceLattice lattice;
    protected final IndexIteratorSizable indexIterator;
    protected boolean rescalingToFitVolume = true;
    protected final AtomActionTranslateTo atomActionTranslateTo;
    protected final AtomGroupAction atomGroupAction;
    protected MyLattice myLat;
    private int[][] indices = null;
    private AtomAgentManager cellManager;
    private static final long serialVersionUID = 2L;
    

    public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
        sim.getDefaults().atomSize = 5.0;
        Phase phase = new Phase(sim);
        sim.addPhase(phase);
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species);
        int k = 4;
        phase.getAgent(species).setNMolecules(4 * k * k * k);
//        ColorSchemeByType colorScheme = new ColorSchemeByType();
        // CubicLattice lattice = new LatticeCubicBcc();
        BravaisLatticeCrystal lattice = new LatticeCubicFcc();
        // CubicLattice lattice = new LatticeCubicSimple();
        ConfigurationMonoclinicLattice configuration = new ConfigurationMonoclinicLattice(lattice);
        // phase.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
        configuration.initializeCoordinates(phase);
        // etomica.graphics.DisplayPhase display = new
        // etomica.graphics.DisplayPhase(phase);

        SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);
//        ((ColorSchemeByType) ((DisplayPhase) simGraphic.displayList()
//                .getFirst()).getColorScheme()).setColor(species
//                .getMoleculeType(), java.awt.Color.red);
        simGraphic.makeAndDisplayFrame(APP_NAME);
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