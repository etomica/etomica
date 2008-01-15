package etomica.models.water;

import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.ISimulation;
import etomica.species.Species;
import etomica.units.ElectronVolt;

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
    public AtomFactoryWater4P(ISimulation sim, Species species) {
		super(new AtomTypeMolecule(species, new AtomPositionGeometricCenter(sim.getSpace())));
        
		AtomWater4P.Echarge[AtomWater4P.indexH1] = ElectronVolt.UNIT.toSim( 0.52);
		AtomWater4P.Echarge[AtomWater4P.indexH2] = ElectronVolt.UNIT.toSim( 0.52);
		AtomWater4P.Echarge[AtomWater4P.indexO] = ElectronVolt.UNIT.toSim( 0.00);
		AtomWater4P.Echarge[AtomWater4P.indexM] = ElectronVolt.UNIT.toSim(-1.04);
		
        AtomTypeSphere hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
        AtomTypeSphere oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.154);
        AtomTypeSphere mType = new AtomTypeSphere(new ElementSimple("M", 1.0), 2.0);
        ((AtomTypeMolecule)atomType).addChildType(hType);
        ((AtomTypeMolecule)atomType).addChildType(oType);
        ((AtomTypeMolecule)atomType).addChildType(mType);

        hFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), hType) :
                                     new AtomFactoryMono(sim.getSpace(), hType);
		oFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), oType) :
		                             new AtomFactoryMono(sim.getSpace(), oType);
		mFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), mType) :
		                             new AtomFactoryMono(sim.getSpace(), mType);
		
		((AtomTypeMolecule)atomType).setConformation(new ConformationWaterTIP4P(sim.getSpace())); 
	}

	/**
	 * @see etomica.atom.AtomFactory#build(etomica.Atom)
	 */
	public IAtom makeAtom() {
        isMutable = false;
        AtomWater4P water = new AtomWater4P(atomType);
		
        water.H1 = (IAtomPositioned)hFactory.makeAtom();
        water.H2 = (IAtomPositioned)hFactory.makeAtom();
        water.O = (IAtomPositioned)oFactory.makeAtom();
        water.M = (IAtomPositioned)mFactory.makeAtom();
        
        water.addChildAtom(water.H1);
        water.addChildAtom(water.H2);
        water.addChildAtom(water.O);
        water.addChildAtom(water.M);
        ((AtomTypeMolecule)atomType).getConformation().initializePositions(water.getChildList());
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
