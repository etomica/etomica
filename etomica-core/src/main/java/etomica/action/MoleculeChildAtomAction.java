/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import java.io.Serializable;

import etomica.atom.IAtomList;
import etomica.atom.IMolecule;

/**
 * Wraps an AtomAction, and performs the wrapped action on the atom
 * only if it is a leaf atom; if given an atom group (as indicated
 * by the atom's node), performs action instead on all the atom's
 * child atoms.  This process continues recursively until the leaf atoms are
 * encountered.
 *
 * @author David Kofke
 */
public class MoleculeChildAtomAction implements MoleculeAction, Serializable {

    /**
     * Constructor takes wrapped action, which is final.
     */
    public MoleculeChildAtomAction(AtomAction action) {
        this.action = action;
    }
    /* (non-Javadoc)
     * @see etomica.action.AtomAction#actionPerformed(etomica.Atom)
     */
    public void actionPerformed(IMolecule atom) {
        IAtomList atomList = atom.getChildList();
        int size = atomList.getAtomCount();
        for(int i=0; i<size; i++) {
            action.actionPerformed(atomList.getAtom(i));
        }
    }

    /**
     * @return Returns the wrapped action.
     */
    public AtomAction getAtomAction() {
        return action;
    }

    private static final long serialVersionUID = 1L;
    private final AtomAction action;
}
