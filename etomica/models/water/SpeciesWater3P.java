package etomica.models.water;
import java.lang.reflect.Constructor;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IMolecule;
import etomica.atom.Molecule;
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
        this(space, new AtomTypeMolecule(new AtomPositionGeometricCenter(space)));
    }
    
    public SpeciesWater3P(Space space, AtomTypeMolecule moleculeType) {
       super(moleculeType);
       this.space = space;
       hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
       oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.167);
       atomType.addChildType(hType);
       atomType.addChildType(oType);

       atomType.setConformation(new ConformationWater3P(space));
    }
    
    public IMolecule makeMolecule() {
        isMutable = false;
        Molecule water = new Molecule(atomType);
        water.addChildAtom(new AtomLeaf(space, hType));
        water.addChildAtom(new AtomLeaf(space, hType));
        water.addChildAtom(new AtomLeaf(space, oType));
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

    public final static int indexH1 = 0;
    public final static int indexH2 = 1;
    public final static int indexO  = 2;

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomTypeSphere oType, hType;
}