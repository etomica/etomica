package etomica.atom;

import etomica.phase.Phase;
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
    
    public void setPhase(Phase p) {
        list = p.getSpeciesMaster().getLeafList();
    }
    
    /**
     * returns a random atom from the phase's leaf atom list
     */
    public Atom getAtom() {
        if (Debug.ON && list== null) throw new IllegalStateException("must set the phase before calling getAtom");
        return list.get(random.nextInt(list.size()));
    }
    
    private static final long serialVersionUID = 1L;
    protected AtomArrayList list = null;
    protected IRandom random;
}
