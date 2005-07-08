package etomica.atom;

import etomica.Atom;
import etomica.Debug;
import etomica.Phase;

/**
 * AtomSource that returns a completely random molecule.
 */
public class AtomSourceRandomMolecule implements AtomSource, java.io.Serializable {

    public void setPhase(Phase p) {
        phase = p;
    }
    
    /**
     * returns a random molecule from the phase
     */
    public Atom getAtom() {
        if (Debug.ON && phase == null) throw new IllegalStateException("must set the phase before calling getAtom");
        return phase.randomMolecule();
    }
    
    protected Phase phase = null;
}
