package etomica;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListCompound;
import etomica.lattice.Crystal;
import etomica.lattice.IndexIteratorSequential;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.LatticeCrystal;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.Primitive;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.space.Vector;
import etomica.space3d.Space3D;

/**
 * Creates a configuration using a CubicLattice to specify positions.  Has
 * capability to assign lattice site to atoms when specifying their coordinates.
 * See setAssigningSitesToAtoms method.
 */
public class ConfigurationLattice extends Configuration implements Atom.AgentSource {

    public ConfigurationLattice(Primitive primitive) {
        this(new LatticeCrystal(new Crystal(primitive, new BasisMonatomic(primitive.space))));
    }
    public ConfigurationLattice(LatticeCrystal lattice) {
        this(lattice, new IndexIteratorSequential(lattice.D()));//need a default iterator
    }
	/**
	 * Constructor for ConfigurationLattice.
	 */
	public ConfigurationLattice(LatticeCrystal lattice, IndexIteratorSizable indexIterator) {
	    super(lattice.getSpace());
        this.lattice = lattice;
        this.indexIterator = indexIterator;
        atomActionTranslateTo = new AtomActionTranslateTo(lattice.getSpace());
        atomIterator = new AtomIteratorListCompound();
        work = space.makeVector();
	}
	
    public void initializePositions(AtomList[] lists) {
        if(lists.length == 0) return;
        atomIterator.setLists(lists);
        int sumOfMolecules = atomIterator.size();
        if(sumOfMolecules == 0) {return;}
        if(sumOfMolecules == 1) {
            atomIterator.reset();
            work.E(0.0);
            atomActionTranslateTo.setDestination(work);
            atomActionTranslateTo.actionPerformed(atomIterator.nextAtom());
            return;
        }
        int basisSize = lattice.getCrystal().getBasis().size();
        int nCells = (int)Math.ceil((double)sumOfMolecules/(double)basisSize);
        
        //determine scaled shape of simulation volume
        double[] shape = (double[])dimensions.clone();
        double[] latticeConstant = lattice.getCrystal().getLattice().getPrimitive().getSize();
        for(int i=0; i<shape.length; i++) {
            shape[i] /= dimensions[0]*latticeConstant[i];
        }
        
        //determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        int[] iteratorDimensions = new int[latticeDimensions.length + 1];
        System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0, latticeDimensions.length);
        iteratorDimensions[latticeDimensions.length] = basisSize;
        indexIterator.setSize(iteratorDimensions);
    
        //determine lattice constant
        Vector latticeScaling = space.makeVector();
        if(rescalingToFitVolume) {
            for(int i=0; i<latticeDimensions.length; i++) {
                //in favorable situations, this should be approximately equal to 1.0
                latticeScaling.setX(i,dimensions[i]/(latticeDimensions[i]*latticeConstant[i]));
            }
        }
        else {
            latticeScaling.E(1.0);
        }

        //determine amount to shift lattice so it is centered in volume
        double[] offset = (double[])dimensions.clone();
        double[] vectorOfMax = new double[lattice.getSpace().D()]; 
        double[] vectorOfMin = new double[lattice.getSpace().D()]; 
        for(int i=0; i<lattice.getSpace().D(); i++) {
            vectorOfMax[i] = Double.NEGATIVE_INFINITY;
            vectorOfMin[i] = Double.POSITIVE_INFINITY;
        }

        //XXX this looks scary and asking for trouble
        indexIterator.reset();
        while(indexIterator.hasNext()) {
            Vector site = (Vector)lattice.site(indexIterator.next());
            site.TE(latticeScaling);
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
        atomIterator.reset();
        indexIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            if (!a.node.isLeaf()) {
                //initialize coordinates of child atoms
                Conformation config = a.type.creator().getConformation();
                config.initializePositions(((AtomTreeNodeGroup)a.node).childList);
            }
            Vector site = (Vector)lattice.site(indexIterator.next());
            site.TE(latticeScaling);
            site.PE(offsetVector);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
            if(assigningSitesToAtoms) ((Agent)a.allatomAgents[siteIndex]).site = (Vector)site.clone();//assign site to atom if so indicated
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
            for (int idim=0; idim<shape.length; idim++) {
                if (latticeDimensions[idim] > 0) continue;
                if (shape[idim] < smin) {
                    smin = shape[idim];
                    dmin = idim;
                }
                product *= shape[idim];
            }
            // round off except for last dimension (then round up)
            if (dimLeft > 1) {
                latticeDimensions[dmin] = (int)Math.round(shape[dmin]*Math.pow((nCellsLeft/product),1.0/dimLeft));
            }
            else {
                latticeDimensions[dmin] = (int)Math.ceil(shape[dmin]*nCellsLeft/product);
            }
            nCellsLeft = (nCellsLeft + latticeDimensions[dmin] - 1) / latticeDimensions[dmin];
            dimLeft--;
        }
        return latticeDimensions;
    }
	
	private final LatticeCrystal lattice;
    private final IndexIteratorSizable indexIterator;
    private final Vector work;
	private boolean rescalingToFitVolume = true;
	private boolean assigningSitesToAtoms = false;
	private int siteIndex = -1;
    private final AtomActionTranslateTo atomActionTranslateTo;
    private final AtomIteratorListCompound atomIterator;

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
        Simulation sim = new Simulation(Space3D.getInstance());
		Default.ATOM_SIZE = 5.0;
		Phase phase = new Phase(sim);
		SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
		int k = 4;
		species.setNMolecules(4*k*k*k);
		etomica.graphics.ColorSchemeByType.setColor(species, java.awt.Color.red);
//        CubicLattice lattice = new LatticeCubicBcc();
        LatticeCrystal lattice = new LatticeCubicFcc();
//        CubicLattice lattice = new LatticeCubicSimple();
		ConfigurationLattice configuration = new ConfigurationLattice(lattice);
//        phase.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
		phase.setConfiguration(configuration);
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
	 * Sets flag that causes sites to be assigned to atoms.  When this is set
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
