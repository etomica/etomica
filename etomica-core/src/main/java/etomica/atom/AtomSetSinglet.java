/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


/**
 * Data structure that contains a single mutable atom instance.
 */
public class AtomSetSinglet implements IAtomList, java.io.Serializable {

    public AtomSetSinglet() {
    }
    
    public AtomSetSinglet(IAtom atom) {
        this.atom = atom;
    }
    
    public final IAtom getAtom(int i) {
        if(i == 0) return atom;
        throw new IllegalArgumentException();
    }

    public final int getAtomCount() {
        return 1;
    }
    
    public String toString() {
        return "["+atom+"]";
    }

    private static final long serialVersionUID = 1L;
    public IAtom atom;
}
