package etomica.models.water;

import etomica.atom.Atom;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactory;
import etomica.space.CoordinateFactorySphere;

/**
 * Factory that constructs a 4-point water molecule, with three child atoms of 
 * two Hydrogen, one Oxygen and an additional point charge.  TIP4P conformation
 * is used by default, but can be changed via setConfiguration.
 * @author kofke
 *
 */
public class AtomFactoryWater4P extends AtomFactory {

	/**
	 * Constructor for AtomFactoryWater.
	 * @param sim
	 */
    public AtomFactoryWater4P(Simulation sim, AtomTypeGroup agentType) {
		super(new AtomTypeGroup(new AtomPositionGeometricCenter(sim.space)), AtomTreeNodeWater4P.FACTORY);
		atomType.setParentType(agentType);
        
        AtomTypeSphere hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
        AtomTypeSphere oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.154);
        AtomTypeSphere mType = new AtomTypeSphere(new ElementSimple("M", 1.0), 2.0);
        hType.setParentType((AtomTypeGroup)atomType);
        oType.setParentType((AtomTypeGroup)atomType);
        mType.setParentType((AtomTypeGroup)atomType);

        CoordinateFactory leafCoordFactory = new CoordinateFactorySphere(sim);
        hFactory = new AtomFactoryMono(leafCoordFactory, hType);
		oFactory = new AtomFactoryMono(leafCoordFactory, oType);
		mFactory = new AtomFactoryMono(leafCoordFactory, mType);
		
		conformation = new ConformationWaterTIP4P(sim.space); 
	}

	/**
	 * @see etomica.atom.AtomFactory#build(etomica.Atom)
	 */
	public Atom makeAtom() {
        Atom group = newParentAtom();
		AtomTreeNodeWater4P waterNode = (AtomTreeNodeWater4P)group.node;
		waterNode.O = (AtomLeaf)oFactory.makeAtom();
        waterNode.H1 = (AtomLeaf)hFactory.makeAtom();
        waterNode.H2 = (AtomLeaf)hFactory.makeAtom();
        waterNode.M = (AtomLeaf)mFactory.makeAtom();
        waterNode.O.node.setParent(waterNode);
        waterNode.H1.node.setParent(waterNode);
        waterNode.H2.node.setParent(waterNode);
        waterNode.M.node.setParent(waterNode);
		conformation.initializePositions(waterNode.childList);
		return group;
	}
    
    /**
     * Returns 5, equal to 1 parent molecule + 4 child atoms in the molecule.
     */
    public int getNumTreeAtoms() {
        return 5;
    }
    
    /**
     * Returns 4.
     */
    public int getNumChildAtoms() {
        return 4;
    }
    
    /**
     * Returns 4.
     */
    public int getNumLeafAtoms() {
        return 4;
    }
    
	public final AtomFactoryMono hFactory, oFactory, mFactory;
}
