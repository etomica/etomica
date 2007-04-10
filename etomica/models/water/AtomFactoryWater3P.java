package etomica.models.water;

import etomica.atom.Atom;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.Simulation;

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
    public AtomFactoryWater3P(Simulation sim) {
		super(new AtomTypeGroup(new AtomPositionGeometricCenter(sim.getSpace())));

        AtomTypeSphere hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
        AtomTypeSphere oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.167);
        hType.setParentType((AtomTypeGroup)atomType);
        oType.setParentType((AtomTypeGroup)atomType);
        hFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), hType) :
                                     new AtomFactoryMono(sim.getSpace(), hType);
		oFactory = sim.isDynamic() ? new AtomFactoryMonoDynamic(sim.getSpace(), oType) :
                                     new AtomFactoryMono(sim.getSpace(), oType);

		conformation = new ConformationWater3P(sim.getSpace()); 
	}

	/**
	 * @see etomica.atom.AtomFactory#build(etomica.Atom)
	 */
	public Atom makeAtom() {
        isMutable = false;
        AtomWater3P water = new AtomWater3P(atomType);
		water.O = (AtomLeaf)oFactory.makeAtom();
        water.H1 = (AtomLeaf)hFactory.makeAtom();
        water.H2 = (AtomLeaf)hFactory.makeAtom();
        water.O.setParent(water);
        water.H1.setParent(water);
        water.H2.setParent(water);
		conformation.initializePositions(water.getChildList());
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
