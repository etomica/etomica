package etomica;

import etomica.lattice.*;

/**
 * @author kofke
 *
 * Creates a configuration using a BravaisLattice to specify positions.  Has
 * capability to assign lattice site to atoms when specifying their coordinates.
 * See setAssigningSitesToAtoms method.
 */
public class ConfigurationLattice extends Configuration implements Atom.AgentSource {

	/**
	 * Constructor for ConfigurationLattice.
	 * @param space
	 */
	public ConfigurationLattice(Space space, BravaisLattice lattice) {
		super(space);
		this.lattice = lattice;
	}
	
	public ConfigurationLattice(Space space, Crystal crystal) {
		this(space, BravaisLattice.makeLattice(space, crystal));
	}

	/**
	 * Constructor for ConfigurationLattice.
	 * @param sim
	 */
	public ConfigurationLattice(Simulation sim, BravaisLattice lattice) {
		super(sim);
		this.lattice = lattice;
	}

	/**
	 * @see etomica.Configuration#initializePositions(etomica.AtomIterator)
	 */
	public void initializePositions(AtomIterator[] iterators) {
		if(iterators == null || iterators.length == 0) return;
		AtomIterator iterator;
		if(iterators.length == 1) iterator = iterators[0];
		else iterator = new AtomIteratorCompound(iterators);//lump 'em all together
		if(iterator.size() == 0) {return;}
		if(iterator.size() == 1) {
			iterator.reset();
			iterator.next().coord.translateTo(space.origin());
			return;
		}
    
	// Count number of molecules
		int sumOfMolecules = iterator.size();
		if(sumOfMolecules == 0) {return;}
		
		if(rescalingToFitVolume) lattice.setSize(sumOfMolecules, dimensions);
		else lattice.setSize(sumOfMolecules);
        
        AtomIteratorList siteIterator = new AtomIteratorList(lattice.siteList());
        
		Space.Vector center = Space.makeVector(dimensions);
		center.TE(0.5);
		lattice.coord.translateCOMTo(center);
		lattice.eventManager.fireEvent((LatticeEvent)null);

   // Place molecules  
		iterator.reset();
		while(iterator.hasNext()) {
			Atom a = iterator.next();
			//initialize coordinates of child atoms
			try {//may get null pointer exception when beginning simulation
				a.creator().getConfiguration().initializeCoordinates(a);
			} catch(NullPointerException e) {}
			Atom site = siteIterator.next();
			a.coord.translateTo(site.coord.position());//use translateTo instead of E because atom might be a group
			if(assigningSitesToAtoms) ((Agent)a.allatomAgents[siteIndex]).site = site;//assign site to atom if so indicated
		}
	}
	
	private BravaisLattice lattice;
	private boolean rescalingToFitVolume = true;
	private boolean assigningSitesToAtoms = false;
	private int siteIndex = -1;

	public static void main(String[] args) {
		etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new Space3D());
		Default.ATOM_SIZE = 5.0;
		Space space = sim.space;
		Phase phase = new Phase(sim);
		SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
		int k = 2;
		species.setNMolecules(4*k*k*k);
		etomica.graphics.ColorSchemeByType.setColor(species, java.awt.Color.red);
		Crystal crystal = new CrystalFcc(space);
		ConfigurationLattice configuration = new ConfigurationLattice(space, crystal);
		phase.setConfiguration(configuration);
		etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase(sim);
		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
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
		public Atom site;
	}

}
