package etomica.models.water;

import etomica.atom.Atom;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomLinker;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactory;
import etomica.space.CoordinateFactoryNull;
import etomica.space.CoordinateFactorySphere;
import etomica.species.Species;

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
		super(new CoordinateFactoryNull(), new AtomTypeGroup(agentType,new AtomPositionGeometricCenter(sim.space)), AtomTreeNodeWater4P.FACTORY);

        AtomTypeSphere hType = new AtomTypeSphere((AtomTypeGroup)atomType, 1.0, 2.0);
        AtomTypeSphere oType = new AtomTypeSphere((AtomTypeGroup)atomType, 16.0, 3.154);
        AtomTypeSphere mType = new AtomTypeSphere((AtomTypeGroup)atomType, 1.0, 2.0);

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
		waterNode.O = oFactory.makeAtom();
        waterNode.H1 = hFactory.makeAtom();
        waterNode.H2 = hFactory.makeAtom();
        waterNode.M = mFactory.makeAtom();
        waterNode.O.node.setParent(waterNode);
        waterNode.H1.node.setParent(waterNode);
        waterNode.H2.node.setParent(waterNode);
        waterNode.M.node.setParent(waterNode);
		conformation.initializePositions(waterNode.childList);
		return group;
	}
    
    public void setSpecies(Species species) {
        atomType.setSpecies(species);
        hFactory.setSpecies(species);
        oFactory.setSpecies(species);
        mFactory.setSpecies(species);
    }

	public final AtomFactoryMono hFactory, oFactory, mFactory;
}
