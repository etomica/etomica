/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


/**
 * Data structure that contains two mutable atom instances.
 */
public class MoleculePair implements IMoleculeList, java.io.Serializable {

    public MoleculePair() {
    }
    
    public MoleculePair(IMolecule atom0, IMolecule atom1) {
        this.atom0 = atom0;
        this.atom1 = atom1;
    }
    
    /* (non-Javadoc)
     * @see etomica.AtomSet#getAtom(int)
     */
    public final IMolecule getMolecule(int i) {
        if(i == 0) return atom0;
        if(i == 1) return atom1;
        throw new IllegalArgumentException();
    }

    /* (non-Javadoc)
     * @see etomica.AtomSet#count()
     */
    public final int getMoleculeCount() {
        return 2;
    }
    
    public String toString() {
        return "["+atom0+","+atom1+"]";
    }

    private static final long serialVersionUID = 1L;
    public IMolecule atom0, atom1;
}
