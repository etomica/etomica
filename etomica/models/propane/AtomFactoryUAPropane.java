package etomica.models.propane;

import etomica.atom.Atom;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactory;
import etomica.space.CoordinateFactorySphere;
import etomica.species.Species;

/**
 * Factory that constructs a 3-point water molecule, with three child atoms of 
 * two Hydrogen and one Oxygen.
 * @author kofke
 *
 */
public class AtomFactoryUAPropane extends AtomFactory {

	/**
	 * Constructor for AtomFactoryWater.
	 * @param sim
	 * @param sequencerFactory
	 */
    public AtomFactoryUAPropane(Simulation sim, AtomTypeGroup agentType) {
		super(new AtomTypeGroup(new AtomPositionGeometricCenter(sim.space)), AtomTreeNodeUAPropane.FACTORY);
        atomType.setParentType(agentType);
        AtomTypeSphere UAType = new AtomTypeSphere(new ElementSimple("UA", 15), 3.75);
        UAType.setParentType((AtomTypeGroup)atomType);
//        AtomTypeSphere oType = new AtomTypeSphere((AtomTypeGroup)atomType, 16.0, 3.167);
        CoordinateFactory leafCoordFactory = new CoordinateFactorySphere(sim);
        UAFactory = new AtomFactoryMono(leafCoordFactory, UAType);
//		oFactory = new AtomFactoryMono(leafCoordFactory, oType);

		conformation = new ConformationUAPropane(sim.space); 
	}

	/**
	 * @see etomica.atom.AtomFactory#build(etomica.Atom)
	 */
	public Atom makeAtom() {
        Atom group = newParentAtom();
		AtomTreeNodeUAPropane propaneNode = (AtomTreeNodeUAPropane)group.node;
//		waterNode.O = (AtomLeaf)oFactory.makeAtom();
        propaneNode.UA1 = (AtomLeaf)UAFactory.makeAtom();
        propaneNode.UA2 = (AtomLeaf)UAFactory.makeAtom();
//        waterNode.O.node.setParent(waterNode);
        propaneNode.UA1.node.setParent(propaneNode);
        propaneNode.UA2.node.setParent(propaneNode);
		conformation.initializePositions(propaneNode.childList);
		return group;
	}
    
    public void setSpecies(Species species) {
        atomType.setSpecies(species);
//        oFactory.setSpecies(species);
    }
    
	public final AtomFactoryMono UAFactory;

	/* (non-Javadoc)
	 * @see etomica.atom.AtomFactory#getNumTreeAtoms()
	 */
	public int getNumTreeAtoms() {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see etomica.atom.AtomFactory#getNumChildAtoms()
	 */
	public int getNumChildAtoms() {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see etomica.atom.AtomFactory#getNumLeafAtoms()
	 */
	public int getNumLeafAtoms() {
		// TODO Auto-generated method stub
		return 0;
	}
}
