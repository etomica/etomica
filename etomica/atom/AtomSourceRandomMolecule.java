package etomica.atom;

import etomica.box.Box;
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
    
    public void setBox(Box p) {
        box = p;
    }
    
    /**
     * returns a random molecule from the box
     */
    public IAtom getAtom() {
        return box.molecule(random.nextInt(box.moleculeCount()));
    }
    
    private static final long serialVersionUID = 1L;
    protected Box box = null;
    protected IRandom random;
}
