package etomica.models.water;

import etomica.api.IMolecule;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeSphere;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species for 4-point water molecule.
 */
public class SpeciesWater4P extends Species {

    public SpeciesWater4P(Space space) {
        super(new AtomPositionGeometricCenter(space));
        this.space = space;
        hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
        oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.167);
        mType = new AtomTypeSphere(new ElementSimple("M", 1.0), 2.0);
        addChildType(hType);
        addChildType(oType);
        addChildType(mType);

        setConformation(new ConformationWaterTIP4P(space)); 
     }

     public IMolecule makeMolecule() {
         Molecule water = new Molecule(this);
         water.addChildAtom(new AtomLeaf(space, hType));
         water.addChildAtom(new AtomLeaf(space, hType));
         water.addChildAtom(new AtomLeaf(space, oType));
         water.addChildAtom(new AtomLeaf(space, mType));
         conformation.initializePositions(water.getChildList());
         return water;
     }

     public AtomTypeSphere getHydrogenType() {
         return hType;
     }

     public AtomTypeSphere getOxygenType() {
         return oType;
     }

     public AtomTypeSphere getMType() {
         return mType;
     }

     public int getNumLeafAtoms() {
         return 4;
     }
    
    public final static int indexH1 = 0;
    public final static int indexH2 = 1;
    public final static int indexO  = 2;
    public final static int indexM  = 3;

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomTypeSphere oType, hType, mType;
}