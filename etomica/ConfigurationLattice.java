package etomica;

import etomica.lattice.*;

/**
 * @author kofke
 *
 * Creates a configuration using a BravaisLattice to specify positions.
 */
public class ConfigurationLattice extends Configuration {

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
		
		lattice.setSize(sumOfMolecules, dimensions);
        
        AtomIteratorList siteIterator = new AtomIteratorList(lattice.siteList());
        
		Atom centerAtom = lattice.centerSite();
		centerAtom = centerAtom.node.firstLeafAtom();
		Space.Vector center = Space.makeVector(dimensions);
		center.TE(0.5);
//		center.ME(centerAtom.coord.position());
		lattice.coord.translateCOMTo(center);

   // Place molecules  
		iterator.reset();
		while(iterator.hasNext()) {
			Atom a = iterator.next();
			//initialize coordinates of child atoms
			try {//may get null pointer exception when beginning simulation
				a.creator().getConfiguration().initializeCoordinates(a);
			} catch(NullPointerException e) {}
			a.coord.translateTo(siteIterator.next().coord.position());//use translateTo instead of E because atom might be a group
		}

	}
	
	private BravaisLattice lattice;

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
}
