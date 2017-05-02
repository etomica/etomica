/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * A box event that is somehow related to an atom.  The atom might have been
 * added or removed, or may have a new index.  Details may be determined from
 * the other interfaces implemented by the event object or obtained from
 * calling methods from those interfaces.
 */
public interface IBoxAtomEvent extends IBoxEvent {
    
    /**
     * @return the atom that is related to this event.
     */
    public IAtom getAtom();
}
