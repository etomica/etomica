package etomica;

import etomica.lattice.Crystal;
import etomica.lattice.LatticeEvent;
import etomica.lattice.SpaceLattice;
import etomica.lattice.crystal.BasisMonatomic;

/**
 * @author kofke
 *
 * Creates a configuration using a SpaceLattice to specify positions.  Has
 * capability to assign lattice site to atoms when specifying their coordinates.
 * See setAssigningSitesToAtoms method.
 */
public class ConfigurationLattice extends Configuration implements Atom.AgentSource {

	/**
	 * Constructor for ConfigurationLattice.
	 * @param space
	 */
	public ConfigurationLattice(SpaceLattice lattice) {
		this(new Crystal(lattice, new BasisMonatomic(lattice.getSpace().D())));
	}
	
	public ConfigurationLattice(Crystal crystal) {
        super(crystal.getSpace());
	    this.crystal = crystal;
	}

	/**
	 * @see etomica.Configuration#initializePositions(etomica.AtomIterator)
	 */
	public void initializePositions(AtomIterator[] iterators) {
		if(iterators == null || iterators.length == 0) return;
		AtomIterator iterator;
		if(iterators.length == 1) iterator = iterators[0];
		else iterator = new AtomIteratorCompound(iterators);
        int sumOfMolecules = iterator.size();
        if(sumOfMolecules == 0) {return;}
		if(sumOfMolecules == 1) {
			iterator.reset();
			iterator.nextAtom().coord.translateTo(space.origin());
			return;
		}
    
		if(rescalingToFitVolume) lattice.setSize(sumOfMolecules, dimensions);
		else lattice.setSize(sumOfMolecules);
        
        AtomIteratorList siteIterator = new AtomIteratorList(lattice.siteList());
        
		Space.Vector center = Space.makeVector(dimensions);
		center.TE(0.5);
		lattice.coord.translateCOMTo(center);
		lattice.eventManager.fireEvent((LatticeEvent)null);

   // Place molecules  
		iterator.reset();
        siteIterator.reset();
		while(iterator.hasNext()) {
			Atom a = iterator.nextAtom();
			//initialize coordinates of child atoms
			try {//may get null pointer exception when beginning simulation
				a.creator().getConfiguration().initializeCoordinates(a);
			} catch(NullPointerException e) {}
			Atom site = siteIterator.nextAtom();
			a.coord.translateTo(site.coord.position());//use translateTo instead of E because atom might be a group
			if(assigningSitesToAtoms) ((Agent)a.allatomAgents[siteIndex]).site = site;//assign site to atom if so indicated
		}
	}
	
	private SpaceLattice lattice;
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

    
//	public static void main(String[] args) {
//		etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new Space3D());
//		Default.ATOM_SIZE = 5.0;
//		Space space = sim.space;
//		Phase phase = new Phase(space);
//		SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
//		int k = 2;
//		species.setNMolecules(4*k*k*k);
//		etomica.graphics.ColorSchemeByType.setColor(species, java.awt.Color.red);
//		Crystal crystal = new CrystalFcc(space);
//		ConfigurationLattice configuration = new ConfigurationLattice(space, crystal);
//		phase.setConfiguration(configuration);
//		etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase(phase);
//		sim.elementCoordinator.go();
//		sim.makeAndDisplayFrame();
//	}
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
		public Space.Vector site;
	}

    private Crystal crystal;
}
