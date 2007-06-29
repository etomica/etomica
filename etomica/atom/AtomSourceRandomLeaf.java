package etomica.atom;

import etomica.box.Box;
import etomica.util.Debug;
import etomica.util.IRandom;

/**
 * AtomSource that returns a completely random leaf atom.
 */
public class AtomSourceRandomLeaf implements AtomSource, java.io.Serializable {

    /**
     * Sets the random number generator used to pick atoms
     */
    public void setRandomNumberGenerator(IRandom newRandom) {
        random = newRandom;
    }
    
    /**
     * Returns the random number generator used to pick atoms
     */
    public IRandom getRandomNumberGenerator() {
        return random;
    }
    
    public void setBox(Box p) {
        list = p.getSpeciesMaster().getLeafList();
    }
    
    /**
     * returns a random atom from the box's leaf atom list
     */
    public IAtom getAtom() {
        if (Debug.ON && list== null) throw new IllegalStateException("must set the box before calling getAtom");
        return list.getAtom(random.nextInt(list.getAtomCount()));
    }
    
    private static final long serialVersionUID = 1L;
    protected AtomSet list = null;
    protected IRandom random;
}
