package etomica;

import etomica.atom.iterator.AtomIteratorCompound;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.CubicLattice;
import etomica.lattice.IndexIteratorSequential;
import etomica.lattice.IndexIteratorSizable;
import etomica.space.Vector;

/**
 * Creates a configuration using a CubicLattice to specify positions.  Has
 * capability to assign lattice site to atoms when specifying their coordinates.
 * See setAssigningSitesToAtoms method.
 */
public class ConfigurationLattice extends Configuration implements Atom.AgentSource {

    public ConfigurationLattice(CubicLattice lattice) {
        this(lattice, new IndexIteratorSequential(lattice.D()));//need a default iterator
    }
	/**
	 * Constructor for ConfigurationLattice.
	 * @param space
	 */
	public ConfigurationLattice(CubicLattice lattice, IndexIteratorSizable indexIterator) {
	    this.lattice = lattice;
        this.indexIterator = indexIterator;
	}
	
	/**
	 * @see etomica.Configuration#initializePositions(etomica.AtomIterator)
	 */
	public void initializePositions(AtomIterator[] iterators) {
		if(iterators == null || iterators.length == 0) return;
		AtomIterator iterator = (iterators.length == 1) ?
                                    iterator = iterators[0] :
                                    new AtomIteratorCompound(iterators);
        int sumOfMolecules = iterator.size();
        if(sumOfMolecules == 0) {return;}
		if(sumOfMolecules == 1) {
			iterator.reset();
			iterator.nextAtom().coord.translateTo(lattice.getSpace().origin());
			return;
		}
        int nCells = (int)Math.ceil((double)sumOfMolecules/(double)lattice.getBasisSize());
        
        //determine scaled shape of simulation volume
        double[] shape = (double[])dimensions.clone();
        for(int i=0; i<shape.length; i++) {
            shape[i] /= dimensions[0];
        }
        
        //determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        int[] iteratorDimensions = new int[latticeDimensions.length + 1];
        System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0, latticeDimensions.length);
        iteratorDimensions[latticeDimensions.length] = lattice.getBasisSize();
        indexIterator.setSize(iteratorDimensions);
    
        //determine lattice constant
        double latticeConstant = lattice.getLatticeConstant();
        if(rescalingToFitVolume) {
            double xmin = Double.MAX_VALUE;
            for(int i=0; i<latticeDimensions.length; i++) {
                double x = dimensions[i]/latticeDimensions[i];//best value of lattice constant for this dimension
                xmin = (x < xmin) ? x : xmin;
            }
            latticeConstant = xmin;
        }
        lattice.setLatticeConstant(latticeConstant);

        //determine amount to shift lattice so it is centered in volume
        double[] offset = (double[])dimensions.clone();
//        for(int i=0; i<shape.length; i++) {
//            offset[i] = 0.5*(dimensions[i] - (latticeDimensions[i]-1)*latticeConstant);
//        }
        double[] vectorOfMax = new double[lattice.getSpace().D()]; 
        double[] vectorOfMin = new double[lattice.getSpace().D()]; 
        for(int i=0; i<lattice.getSpace().D(); i++) {
            vectorOfMax[i] = Double.NEGATIVE_INFINITY;
            vectorOfMin[i] = Double.POSITIVE_INFINITY;
        }
        indexIterator.reset();
        while(indexIterator.hasNext()) {
            Vector site = (Vector)lattice.site(indexIterator.next());
            for(int i=0; i<site.D(); i++) {
                vectorOfMax[i] = Math.max(vectorOfMax[i],site.x(i));
                vectorOfMin[i] = Math.min(vectorOfMin[i],site.x(i));
            }
        }
        for(int i=0; i<lattice.getSpace().D(); i++) {
            offset[i] = 0.5*(dimensions[i] - (vectorOfMax[i]-vectorOfMin[i])) - vectorOfMin[i];
        }
        Vector offsetVector = Space.makeVector(offset);
        
        // Place molecules  
		iterator.reset();
        indexIterator.reset();
		while(iterator.hasNext()) {
			Atom a = iterator.nextAtom();
			//initialize coordinates of child atoms
			try {//may get null pointer exception when beginning simulation
				a.creator().getConfiguration().initializeCoordinates(a);
			} catch(NullPointerException e) {}
            Vector site = (Vector)lattice.site(indexIterator.next());
            site.PE(offsetVector);
			a.coord.translateTo(site);//use translateTo instead of E because atom might be a group
			if(assigningSitesToAtoms) ((Agent)a.allatomAgents[siteIndex]).site = site;//assign site to atom if so indicated
		}
	}
    
    public double getLatticeConstant() {
        return lattice.getLatticeConstant();
    }
    
    private int[] calculateLatticeDimensions(int nCells, double[] shape) {
        double product = 1.0;
        for(int i=0; i<shape.length; i++) product *= shape[i];
        int n = (int)Math.ceil(Math.pow(nCells/product,1.0/shape.length));
        int[] latticeDimensions = new int[shape.length];
        for(int i=0; i<shape.length; i++) {
            latticeDimensions[i] = (int)Math.round((n*shape[i]));//(int)Math.ceil(n*shape[i]);
        }
        return latticeDimensions;
    }
	
	private CubicLattice lattice;
    private IndexIteratorSizable indexIterator;
	private boolean rescalingToFitVolume = true;
	private boolean assigningSitesToAtoms = false;
	private int siteIndex = -1;

//  /**
//  * Sets the size of the lattice (number of atoms in each direction) so that
//  * it has a number of sites equal or greater than the given value n. Sets
//  * lattice to be square (same number in all directions), and finds smallest
//  * size that gives number of sites equal or exceeding the given value.
//  */
// public void setSize(int n) {
//     if (n < 0 || n >= Integer.MAX_VALUE)
//         throw new IllegalArgumentException(
//                 "Inappropriate size specified for lattice: " + n);
//     int i = (int)Math.pow(n,1/D());
//     while(i <= n) {
//         if ((int)Math.pow(i,D()) >= n) break;
//         i++;
//     }
//     setSize(primitive.space.makeArrayD(i));
// }
//
// /**
//  * Performs same actions as setSize(int), then size of primitive vectors are
//  * adjusted such that lattice will fit in the given dimensions. Assumes all
//  * dimension values are equal (future development will address case of non-
//  * square dimensions).
//  */
// public void setSize(int n, double[] dimensions) {
//     if (n < 0 || n >= Integer.MAX_VALUE)
//         throw new IllegalArgumentException(
//                 "Inappropriate size specified for lattice: " + n);
//     int i = (int)Math.pow(n,1/D());
//     while(i <= n) {
//         if ((int)Math.pow(i,D()) >= n) break;
//         i++;
//     }
//     //size primitive vectors to given dimensions
//     double[] newLength = new double[D()];
//     for (int j = 0; j < D(); j++)
//         newLength[j] = dimensions[j] / (double) i;
//     primitive.setSize(newLength);
//     setSize(primitive.space.makeArrayD(i));
// }

// /**
//  * Translates the lattice so that the first site is positioned at the
//  * given point.
//  */
// public void shiftFirstTo(Space.Vector r) {
//     Space.Vector[] coords = (Space.Vector[])sites();
//     for(int i=coords.length-1; i>=0; i--) {
//         coords[i].ME(r);
//     }
// }

    
	public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.INSTANCE);
		Default.ATOM_SIZE = 5.0;
		Space space = sim.space;
		Phase phase = new Phase(space);
		SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
		int k = 4;
		species.setNMolecules(4*k*k*k);
		etomica.graphics.ColorSchemeByType.setColor(species, java.awt.Color.red);
//        CubicLattice lattice = new LatticeCubicBcc();
        CubicLattice lattice = new LatticeCubicFcc();
//        CubicLattice lattice = new LatticeCubicSimple();
		ConfigurationLattice configuration = new ConfigurationLattice(lattice);
//        phase.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
		phase.setConfiguration(configuration);
        phase.speciesMaster.addSpecies(species);
//		etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase(phase);
		
        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(sim);
		simGraphic.makeAndDisplayFrame();
	}
	/**
	 * Returns the resizeLatticeToFitVolume flag.
	 * @return boolean
	 */
	public boolean isRescalingToFitVolume() {
		return rescalingToFitVolume;
	}

	/**
	 * Sets the resizeLatticeToFitVolume flag, which if true indicates that the
	 * primitive vectors should be resized to fit the dimensions of the phase.
	 * Default is true.
	 * @param resizeLatticeToFitVolume The resizeLatticeToFitVolume to set
	 */
	public void setRescalingToFitVolume(boolean resizeLatticeToFitVolume) {
		this.rescalingToFitVolume = resizeLatticeToFitVolume;
	}

	/**
	 * Returns the assigningSitesToAtoms field.
	 * @return boolean
	 */
	public boolean isAssigningSitesToAtoms() {
		return assigningSitesToAtoms;
	}

	/**
	 * Sets flage that causes sites to be assigned to atoms.  When this is set
	 * to true, an atom agent is made in every new Atom instance, which will
	 * hold a field that points to the site used to set the atom's position
	 * during any call to initializeCoordinates.
	 * @param assigningSitesToAtoms The assigningSitesToAtoms to set
	 */
	public void setAssigningSitesToAtoms(boolean assigningSitesToAtoms) {
		this.assigningSitesToAtoms = assigningSitesToAtoms;
		if(assigningSitesToAtoms && siteIndex < 0) siteIndex = Atom.requestAgentIndex(this);
	}
	
	/**
	 * Returns the index used to access the site.  Given an atom, its site is
	 * accessed via: 
	 * ((ConfigurationLattice.Agent)atom.allAtomAgents[siteIndex]).site
	 * @return int the site index in allatomAgents with this class' agent
	 */
	public final int siteIndex() {return siteIndex;}
	
	/**
	 * Implementation of Atom.AgentSource interface.
	 * @see etomica.Atom.AgentSource#makeAgent(Atom)
	 */
	public Object makeAgent(Atom a) {return new Agent();}
	
	/**
	 * Atom agent that simply has a field that points to a lattice site
	 * associated with the atom.
	 */
	public static class Agent {
		public Vector site;
	}

}
