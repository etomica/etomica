package etomica.models.water;

import etomica.Atom;
import etomica.Space;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;

/**
 * Factory that constructs a water molecule, with three child atoms of two
 * Hydrogen and one Oxygen.
 * @author kofke
 *
 */
public class AtomFactoryWater extends AtomFactory {

	/**
	 * Constructor for AtomFactoryWater.
	 * @param sim
	 * @param sequencerFactory
	 */
	public AtomFactoryWater(Space space) {
        this(space, AtomSequencerFactory.SIMPLE);
    }
    
    public AtomFactoryWater(Space space, AtomSequencerFactory sequencerFactory) {
		super(space, sequencerFactory, AtomTreeNodeWater.FACTORY);

		hFactory = new AtomFactoryMono(space, AtomSequencerFactory.SIMPLE);
		oFactory = new AtomFactoryMono(space, AtomSequencerFactory.SIMPLE);
		AtomType hType = new AtomTypeSphere(hFactory, 1.0, /*Electron.UNIT.toSim(0.4238),*/ 2.0);
		AtomType oType = new AtomTypeSphere(oFactory, 16.0, /*Electron.UNIT.toSim(-0.8476),*/ 3.167);
        
		hFactory.setType(hType);
		oFactory.setType(oType);

		configuration = new ConfigurationWater(space); 
	}

	/**
	 * @see etomica.atom.AtomFactory#build(etomica.Atom)
	 */
	public Atom makeAtom() {
        Atom group = newParentAtom();
		AtomTreeNodeWater waterNode = (AtomTreeNodeWater)group.node;
		waterNode.O = oFactory.makeAtom();
        waterNode.H1 = hFactory.makeAtom();
        waterNode.H2 = hFactory.makeAtom();
        waterNode.O.node.setParent(waterNode);
        waterNode.H1.node.setParent(waterNode);
        waterNode.H2.node.setParent(waterNode);
		configuration.initializeCoordinates(group);
		return group;
	}

	/**
	 * @see etomica.atom.AtomFactory#isGroupFactory()
	 */
	public boolean isGroupFactory() {
		return true;
	}

	public final AtomFactoryMono hFactory, oFactory;
}
