/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.util.Debug;

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
        list = p.getLeafList();
//        System.out.println("setBox for getAtom"+p);
    }
    
    /**
     * returns a random atom from the box's leaf atom list
     */
    public IAtom getAtom() {
        if (Debug.ON && list== null) throw new IllegalStateException("must set the box before calling getAtom");
        int n = list.getAtomCount();
        if (n == 0) return null;
        IAtom returnAtom = list.getAtom(random.nextInt(n));
        return returnAtom;
    }
    
    private static final long serialVersionUID = 1L;
    protected IAtomList list = null;
    protected IRandom random;
}
