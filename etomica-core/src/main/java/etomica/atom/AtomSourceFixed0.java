package etomica.atom;

import etomica.box.Box;
import etomica.util.Debug;
import etomica.util.random.IRandom;

public class AtomSourceFixed0 implements AtomSource {

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

    @Override
    public void setBox(Box p) {
        list = p.getLeafList();
    }

    @Override
    public IAtom getAtom() {
        if (Debug.ON && list == null) throw new IllegalStateException("must set the box before calling getAtom");
        int n = list.size();
        if (n == 0) return null;
        return list.get(1+random.nextInt(n-1)); // Fix atom0 by skipping it (Einstein Molecule!) //1,2,...,n-1
    }

    protected IAtomList list = null;
    protected IRandom random;
}
