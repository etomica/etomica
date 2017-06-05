/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

/**
 * A box listener receives notifications for box events.  Each type of event
 * triggers a call to a different method.  A listener receives event
 * notifications once it is added to the box's event manager.
 */
public interface BoxEventListener {

    /**
     * Called when a molecule is added to the box.
     *
     * @param e the event, which contains the added molecule and the box to
     *          which it was added.
     */
    void boxMoleculeAdded(BoxMoleculeEvent e);

    /**
     * Called when a molecule is removed from the box.
     *
     * @param e the event, which contains the removed molecule and the box from
     *          which it was removed.
     */
    void boxMoleculeRemoved(BoxMoleculeEvent e);

    /**
     * Called when the maximum leaf index in the box has changed.
     *
     * @param e the event, which contains the box and the new maximum global
     *          index
     */
    void boxGlobalAtomLeafIndexChanged(BoxIndexEvent e);

    /**
     * Called when an atom's global (leaf) index has changed.  This typically
     * happens after one molecule is removed and the indices of other atoms are
     * renumbered to fill in the gap left by the removed atoms.
     *
     * @param e the event, which contains the box, the atom whose index changed
     *          and the atom's old leaf index.  The atom's new leaf index can
     *          be obtained from the atom.
     */
    void boxAtomLeafIndexChanged(BoxAtomIndexEvent e);

    /**
     * Called when an molecule's index has changed.  This typically happens
     * after one molecule is removed and the indices of another molecule is
     * changed to fill in the gap.
     *
     * @param e the event, which contains the box, the molecule whose index
     *          changed and the molecule's old index.  The molecule's new index
     *          can be obtained from the atom.
     */
    void boxMoleculeIndexChanged(BoxMoleculeIndexEvent e);

    /**
     * Called when the number of molecules of a particular species has changed
     * in the box.  This method is called before the new molecules have
     * been added to the box.
     *
     * @param e the event, which contains the box, the species and the new
     *          number of molecules of that species.
     */
    void boxNumberMolecules(BoxMoleculeCountEvent e);

}
