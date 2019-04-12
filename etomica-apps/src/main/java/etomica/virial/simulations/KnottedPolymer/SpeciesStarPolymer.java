package etomica.virial.simulations.KnottedPolymer;

import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Bead-model species for polymers
 */

public class SpeciesStarPolymer extends Species {

    protected final Space space;
    protected final int nBead = 201; // 201
    protected final AtomType bead;
    protected final AtomType core;
    protected final double mass = 12.0;
    protected final double maxDiameter;
    private boolean isDynamic;


    public SpeciesStarPolymer(Space space, String filename) {
        super();
        this.space = space;
        this.bead = AtomType.simple("Bead", mass);
        this.core = AtomType.simple("Core", mass);
        addChildType(bead);
        addChildType(core);

        ConformationStarPolymer conf = new ConformationStarPolymer(space, filename);
        setConformation(conf);
        this.maxDiameter = conf.maxDiameter;
    }

    @Override
    public IMolecule makeMolecule() {
        Molecule starPolymer = new Molecule(this, nBead);

        if (isDynamic) {
            if (nBead > 1) {
                for (int i = 0; i < nBead - 1; i++) {
                    starPolymer.addChildAtom(new AtomLeafDynamic(space, bead));
                }
                starPolymer.addChildAtom(new AtomLeafDynamic(space, core));
            }
        } else {
            if (nBead > 1) {
                for (int i = 0; i < nBead - 1; i++) {
                    starPolymer.addChildAtom(new Atom(space, bead));
                }
            }
            starPolymer.addChildAtom(new Atom(space, core));
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
}
