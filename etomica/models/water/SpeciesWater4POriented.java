package etomica.models.water;

import etomica.api.IMolecule;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.MoleculeOriented;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.ISpace;
import etomica.species.SpeciesOriented;

public class SpeciesWater4POriented extends SpeciesOriented {

    public SpeciesWater4POriented(ISpace space, boolean isDynamic) {
        super(space);
        this.isDynamic = isDynamic;
        this.space = space;
        hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0);
        oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.167);
        mType = new AtomTypeSphere(Oxygen.INSTANCE, 0);
        addChildType(hType);
        addChildType(oType);
        addChildType(mType);

        setConformation(new ConformationWaterTIP4P(space));
        init();
    }
    
    public IMolecule makeMolecule() {
        MoleculeOriented water = isDynamic ? new MoleculeOrientedDynamic(space, this) :
                                             new MoleculeOriented(space, this);
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
    protected final ISpace space;
    protected final AtomTypeSphere oType, hType, mType;
    protected final boolean isDynamic;
}
