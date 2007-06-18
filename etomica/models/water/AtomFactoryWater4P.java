package etomica.models.water;

import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.ISimulation;

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
    public AtomFactoryWater4P(ISimulation sim) {
		super(new AtomTypeGroup(new AtomPositionGeometricCenter(sim.getSpace())));
        
        AtomTypeSphere hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
        AtomTypeSphere oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.154);
        AtomTypeSphere mType = new AtomTypeSphere(new ElementSimple("M", 1.0), 2.0);
        hType.setParentType((AtomTypeGroup)atomType);
        oType.setParentType((AtomTypeGroup)atomType);
        mType.setParentType((AtomTypeGroup)atomType);

        hFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), hType) :
                                     new AtomFactoryMono(sim.getSpace(), hType);
		oFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), oType) :
		                             new AtomFactoryMono(sim.getSpace(), oType);
		mFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), mType) :
		                             new AtomFactoryMono(sim.getSpace(), mType);
		
		((AtomTypeGroup)atomType).setConformation(new ConformationWaterTIP4P(sim.getSpace())); 
	}

	/**
	 * @see etomica.atom.AtomFactory#build(etomica.Atom)
	 */
	public IAtom makeAtom() {
        isMutable = false;
        AtomWater4P water = new AtomWater4P(atomType);
		water.O = (IAtomPositioned)oFactory.makeAtom();
        water.H1 = (IAtomPositioned)hFactory.makeAtom();
        water.H2 = (IAtomPositioned)hFactory.makeAtom();
        water.M = (IAtomPositioned)mFactory.makeAtom();
        water.addChildAtom(water.O);
        water.addChildAtom(water.H1);
        water.addChildAtom(water.H2);
        water.addChildAtom(water.M);
        ((AtomTypeGroup)atomType).getConformation().initializePositions(water.getChildList());
		return water;
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
    
    private static final long serialVersionUID = 1L;
	public final AtomFactoryMono hFactory, oFactory, mFactory;
}
