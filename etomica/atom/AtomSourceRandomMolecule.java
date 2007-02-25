package etomica.atom;

import etomica.phase.Phase;
import etomica.util.IRandom;

/**
 * AtomSource that returns a completely random molecule.
 */
public class AtomSourceRandomMolecule implements AtomSource, java.io.Serializable {

    /**
     * Sets the random number generator used to pick atoms
     */
    public void setRandom(IRandom newRandom) {
        random = newRandom;
    }
    
    /**
     * Returns the random number generator used to pick atoms
     */
    public IRandom getRandomNumberGenerator() {
        return random;
    }
    
    public void setPhase(Phase p) {
        phase = p;
    }
    
    /**
     * returns a random molecule from the phase
     */
    public Atom getAtom() {
        return phase.molecule(random.nextInt(phase.moleculeCount()));
    }
    
    private static final long serialVersionUID = 1L;
    protected Phase phase = null;
    protected IRandom random;
}
