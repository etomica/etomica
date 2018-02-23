/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;


/**
 * Data structure that contains a single mutable atom instance.
 */
public class MoleculeSetSinglet implements IMoleculeList, java.io.Serializable {

    public MoleculeSetSinglet() {
    }
    
    public MoleculeSetSinglet(IMolecule atom) {
        this.atom = atom;
    }
    
    public final IMolecule get(int i) {
        if(i == 0) return atom;
        throw new IllegalArgumentException();
    }

    public final int size() {
        return 1;
    }
    
    public String toString() {
        return "["+atom+"]";
    }

    private static final long serialVersionUID = 1L;
    public IMolecule atom;
}
