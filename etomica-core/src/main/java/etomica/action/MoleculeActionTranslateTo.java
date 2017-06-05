/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;


import java.io.Serializable;

import etomica.api.IMolecule;
import etomica.atom.IMoleculePositionDefinition;
import etomica.space.Vector;
import etomica.atom.MoleculePositionCOM;
import etomica.space.Space;

/**
 * Moves (translates) an atom to a specified position.  Location of the
 * atom (which may be an atom group) is determined by an AtomPositionDefinition
 * instance that may be set for this class.
 */
public class MoleculeActionTranslateTo implements MoleculeAction, Serializable {
    
    private static final long serialVersionUID = 1L;
    private final Vector destination;
    private IMoleculePositionDefinition atomPositionDefinition;
    private MoleculeChildAtomAction atomTranslator;
    private final Vector translationVector;

    /**
     * Creates new action with atom position defined by its
     * center of mass (via MoleculePositionCOM).
     * @param space
     */
    public MoleculeActionTranslateTo(Space space) {
        destination = space.makeVector();
        atomPositionDefinition = new MoleculePositionCOM(space);
        atomTranslator = new MoleculeChildAtomAction(new AtomActionTranslateBy(space));
        translationVector = ((AtomActionTranslateBy)atomTranslator.getAtomAction()).getTranslationVector();
    }
    
    public void actionPerformed(IMolecule atom) {
        Vector currentPosition = atomPositionDefinition.position(atom);
        translationVector.Ev1Mv2(destination, currentPosition);
        atomTranslator.actionPerformed(atom);
    }
       
    /**
     * @return Returns the destination, the position that the
     * atom will be moved to by this action.
     */
    public Vector getDestination() {
        return destination;
    }
    /**
     * @param newDestination The destination to set.  A local copy
     * is made of the given vector.
     */
    public void setDestination(Vector newDestination) {
        destination.E(newDestination);
    }
    /**
     * @return Returns the atomPositionDefinition.
     */
    public IMoleculePositionDefinition getAtomPositionDefinition() {
        return atomPositionDefinition;
    }
    /**
     * @param atomPositionDefinition The atomPositionDefinition to set.
     */
    public void setAtomPositionDefinition(
            IMoleculePositionDefinition atomPositionDefinition) {
        this.atomPositionDefinition = atomPositionDefinition;
    }
    
    /**
     * Returns the vector that was used to accomplish the most recent translation action.
     * This vector can be used to reverse the translation by multiplying it by -1 and 
     * performing an atomActionTranslateBy with it.
     */
    public Vector getTranslationVector() {
        return translationVector;
    }
}
