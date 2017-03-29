/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import java.io.Serializable;

import etomica.api.IAtom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;

/**
 * 
 * Moves (translates) an atom by a specified vector amount.
 * To move all atoms in a molecule (or atom group), wrap an
 * instance of this class in an AtomGroupAction.
 * 
 * @author David Kofke
 */
public class AtomActionTranslateBy implements AtomAction, Serializable {
    
    private static final long serialVersionUID = 1L;
    private final IVectorMutable translationVector;
    
    public AtomActionTranslateBy(ISpace space) {
        translationVector = space.makeVector();
    }
    
    public void actionPerformed(IAtom atom) {
        atom.getPosition().PE(translationVector);
    }
       
    /**
     * Returns the translation vector, the distance that the
     * atom will be moved by this action. Returns the vector used by this
     * instance, not a copy, so any manipulation of the returned vector will
     * affect the action of this instance.
     */
    public IVectorMutable getTranslationVector() {
        return translationVector;
    }
    /**
     * @param translationVector The translation vector to set.  A local copy
     * is made of the given vector.
     */
    public void setTranslationVector(IVector translationVector) {
        this.translationVector.E(translationVector);
    }
}
