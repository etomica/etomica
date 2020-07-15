package etomica.starpolymer;

import etomica.atom.Atom;
import etomica.atom.AtomType;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Bead-model species for polymers
 */

public class SpeciesStarPolymerSelected extends Species {

    protected final Space space;
    protected final int nBead = 201;
    protected final AtomType bead;
    protected final AtomType core;
    protected final double mass = 12.0;

    public SpeciesStarPolymerSelected(Space space, String filename, int index) {
        super();
        this.space = space;
        this.bead = AtomType.simple("Bead", mass);
        this.core = AtomType.simple("Core", mass);
        addChildType(bead);
        addChildType(core);
        setConformation(new ConformationStarPolymerAll(space, filename, index));
    }

    @Override
    public IMolecule makeMolecule() {
        Molecule starPolymer = new Molecule(this, nBead);

        for (int i = 0; i < nBead - 1; i++) {
            starPolymer.addChildAtom(new Atom(space, bead));
        }
        starPolymer.addChildAtom(new Atom(space, core));

        conformation.initializePositions(starPolymer.getChildList());

        return starPolymer;
    }
}
