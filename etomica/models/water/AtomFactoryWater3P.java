package etomica.models.water;

import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.ISimulation;
import etomica.species.Species;

/**
 * Factory that constructs a 3-point water molecule, with three child atoms of 
 * two Hydrogen and one Oxygen.
 * @author kofke
 *
 */
public class AtomFactoryWater3P extends AtomFactory {

	/**
	 * Constructor for AtomFactoryWater.
	 * @param sim
	 * @param sequencerFactory
	 */
    public AtomFactoryWater3P(ISimulation sim, Species species) {
		super(new AtomTypeMolecule(species, new AtomPositionGeometricCenter(sim.getSpace())));

        AtomTypeSphere hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
        AtomTypeSphere oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.167);
        ((AtomTypeMolecule)atomType).addChildType(hType);
        ((AtomTypeMolecule)atomType).addChildType(oType);
        hFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), hType) :
                                     new AtomFactoryMono(sim.getSpace(), hType);
		oFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), oType) :
                                     new AtomFactoryMono(sim.getSpace(), oType);

		((AtomTypeMolecule)atomType).setConformation(new ConformationWater3P(sim.getSpace())); 
	}

	/**
	 * @see etomica.atom.AtomFactory#build(etomica.Atom)
	 */
	public IAtom makeAtom() {
        isMutable = false;
        AtomWater3P water = new AtomWater3P(atomType);
		water.O = (IAtomPositioned)oFactory.makeAtom();
        water.H1 = (IAtomPositioned)hFactory.makeAtom();
        water.H2 = (IAtomPositioned)hFactory.makeAtom();
        water.addChildAtom(water.O);
        water.addChildAtom(water.H1);
        water.addChildAtom(water.H2);
        ((AtomTypeMolecule)atomType).getConformation().initializePositions(water.getChildList());
		return water;
	}
    
    /**
     * Returns 4, equal to 1 parent molecule + 3 atoms child atoms in the molecule.
     */
    public int getNumTreeAtoms() {
        return 4;
    }
    
    /**
     * Returns 3.
     */
    public int getNumChildAtoms() {
        return 3;
    }
    
    /**
     * Returns 3.
     */
    public int getNumLeafAtoms() {
        return 3;
    }
    
    private static final long serialVersionUID = 1L;
	public final AtomFactoryMono hFactory, oFactory;
}
