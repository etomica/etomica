package etomica.starpolymer;

import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Bead-model species for polymers
 */

public class SpeciesStarPolymerNew extends Species {

    protected final Space space;
    protected final AtomType bead;
    //    protected final AtomType core;
    protected final double mass = 0.00;
    protected int nBead;
    protected boolean isDynamic;
    protected AtomType leafAtomType;
    public int f, l;

    public SpeciesStarPolymerNew(Space space, int f, int l) {
        super();
        this.f = f;
        this.l = l;
        this.space = space;
        this.nBead = f * l + 1;
        this.bead = AtomType.simple("Bead", mass);
//        this.core = AtomType.simple("Core", mass);
        this.leafAtomType = bead;
//        addChildType(core);
        addChildType(bead);

        ConformationStarPolymerGraft conf = new ConformationStarPolymerGraft(space, f, l);
        setConformation(conf);
    }

    @Override
    public IMolecule makeMolecule() {
        Molecule starPolymer = new Molecule(this, nBead);

        if (isDynamic) {
            starPolymer.addChildAtom(new AtomLeafDynamic(space, bead));
        } else {
            starPolymer.addChildAtom(new Atom(space, bead));
        }
        if (nBead > 1) {
            for (int i = 0; i < nBead - 1; i++) {
                starPolymer.addChildAtom(makeLeafAtom());
            }
        }
        conformation.initializePositions(starPolymer.getChildList());
        return starPolymer;
    }

    public void setIsDynamic(boolean newIsDynamic) {
        isDynamic = newIsDynamic;
    }

    public boolean isDynamic() {
        return isDynamic;
    }

    private IAtom makeLeafAtom() {
        return isDynamic ? new AtomLeafDynamic(space, leafAtomType)
                : new Atom(space, leafAtomType);
    }
}
