package etomica.models.water;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeMoleculeOriented;
import etomica.atom.IMolecule;
import etomica.atom.MoleculeOriented;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.space.Space;

public class SpeciesWater3POriented extends SpeciesWater3P {

    public SpeciesWater3POriented(Space space, boolean isDynamic) {
        super(space, new AtomTypeMoleculeOriented(space));
        this.isDynamic = isDynamic;
        ((AtomTypeMoleculeOriented)atomType).init();
    }

    public IMolecule makeMolecule() {
        isMutable = false;
        MoleculeOriented water = isDynamic ? new MoleculeOrientedDynamic(space, atomType) :
                                             new MoleculeOriented(space, atomType);
        water.addChildAtom(new AtomLeaf(space, hType));
        water.addChildAtom(new AtomLeaf(space, hType));
        water.addChildAtom(new AtomLeaf(space, oType));
        atomType.getConformation().initializePositions(water.getChildList());
        return water;
    }

    private static final long serialVersionUID = 1L;
    protected final boolean isDynamic;
}
