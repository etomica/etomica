package etomica.models.water;

import etomica.*;

/**
 * @author kofke
 *
 * Factory that constructs a water molecule, with three child atoms of two
 * Hydrogen and one Oxygen.
 */
public class AtomFactoryWater extends AtomFactory {

	/**
	 * Constructor for AtomFactoryWater.
	 * @param sim
	 * @param sequencerFactory
	 */
	public AtomFactoryWater(Simulation sim) {
		super(sim, sim.getIteratorFactory().neighborSequencerFactory(), AtomTreeNodeWater.FACTORY);

		hFactory = new AtomFactoryMono(simulation, simulation.getIteratorFactory().simpleSequencerFactory());
		oFactory = new AtomFactoryMono(simulation, simulation.getIteratorFactory().simpleSequencerFactory());
		AtomType hType = new AtomType.Sphere(hFactory, 1.0, /*Electron.UNIT.toSim(0.4238),*/ 2.0/2);
		AtomType oType = new AtomType.Sphere(oFactory, 16.0, /*Electron.UNIT.toSim(-0.8476),*/ 3.167/2);
        
		hFactory.setType(hType);
		oFactory.setType(oType);

		configuration = new ConfigurationWater(simulation); 
	}

	/**
	 * @see etomica.AtomFactory#build(etomica.Atom)
	 */
	public Atom build(Atom group) {
		AtomTreeNodeWater waterNode = (AtomTreeNodeWater)group.node;
		waterNode.O = oFactory.makeAtom(waterNode);
		waterNode.H1 = hFactory.makeAtom(waterNode);
		waterNode.H2 = hFactory.makeAtom(waterNode);
		bondInitializer.makeBonds(group);
		configuration.initializeCoordinates(group);
		return group;
	}

	/**
	 * @see etomica.AtomFactory#isGroupFactory()
	 */
	public boolean isGroupFactory() {
		return true;
	}

	private AtomFactoryMono hFactory, oFactory;
}
