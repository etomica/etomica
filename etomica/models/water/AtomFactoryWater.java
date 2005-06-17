package etomica.models.water;

import etomica.Atom;
import etomica.AtomFactory;
import etomica.AtomTypeGroup;
import etomica.Space;
import etomica.Species;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomSequencerFactory;
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
	public AtomFactoryWater(Space space, AtomTypeGroup agentType) {
        this(space, AtomSequencerFactory.SIMPLE, agentType);
    }
    
    public AtomFactoryWater(Space space, AtomSequencerFactory sequencerFactory, AtomTypeGroup agentType) {
		super(space, new AtomTypeGroup(agentType,new AtomPositionGeometricCenter(space)), sequencerFactory, AtomTreeNodeWater.FACTORY);

        AtomTypeSphere hType = new AtomTypeSphere((AtomTypeGroup)atomType, 1.0, /*Electron.UNIT.toSim(0.4238),*/ 2.0);
        AtomTypeSphere oType = new AtomTypeSphere((AtomTypeGroup)atomType, 16.0, /*Electron.UNIT.toSim(-0.8476),*/ 3.167);

        hFactory = new AtomFactoryMono(space, hType, AtomSequencerFactory.SIMPLE);
		oFactory = new AtomFactoryMono(space, oType, AtomSequencerFactory.SIMPLE);

		conformation = new ConformationWater(space); 
	}

	/**
	 * @see etomica.AtomFactory#build(etomica.Atom)
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
		conformation.initializePositions(waterNode.childList);
		return group;
	}
    
    public void setSpecies(Species species) {
        atomType.setSpecies(species);
        hFactory.setSpecies(species);
        oFactory.setSpecies(species);
    }
    
//    public void setDepth(int depth) {
//        atomType.setDepth(depth);
//        hFactory.setDepth(depth+1);
//        oFactory.setDepth(depth+1);
//    }

	public final AtomFactoryMono hFactory, oFactory;
}
