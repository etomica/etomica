/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * Box event that indicates that an atom's index has changed.
 */
public interface IBoxAtomIndexEvent extends IBoxAtomEvent {

    /**
     * @return the index of the atom
     */
    public int getIndex();
}
