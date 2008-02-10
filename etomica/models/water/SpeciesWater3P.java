package etomica.models.water;
import java.lang.reflect.Constructor;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IMolecule;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.ISimulation;
import etomica.space.Space;
import etomica.species.Species;
import etomica.species.SpeciesSignature;

/**
 * Species for 3-point water molecule.
 */
public class SpeciesWater3P extends Species {
    
    public SpeciesWater3P(Space space) {
       super(new AtomTypeMolecule(new AtomPositionGeometricCenter(space)));
       this.space = space;
       hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
       oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.167);
       atomType.addChildType(hType);
       atomType.addChildType(oType);

       atomType.setConformation(new ConformationWater3P(space)); 
    }
    
    public IMolecule makeMolecule() {
        isMutable = false;
        AtomWater3P water = new AtomWater3P(atomType);
        water.H1 = new AtomLeaf(space, hType);
        water.H2 = new AtomLeaf(space, hType);
        water.O = new AtomLeaf(space, oType);
        water.addChildAtom(water.H1);
        water.addChildAtom(water.H2);
        water.addChildAtom(water.O);
        atomType.getConformation().initializePositions(water.getChildList());
        return water;
    }
    
    public AtomTypeSphere getHydrogenType() {
        return hType;
    }
    
    public AtomTypeSphere getOxygenType() {
        return oType;
    }

    public int getNumLeafAtoms() {
        return 3;
    }
    
    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(constructor,new Object[]{});
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomTypeSphere oType, hType;
}