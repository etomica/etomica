/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

/**
 * Box event that indicates the number of molecules of a particular species has
 * changed.  This event is fired before the molecules have actually been added
 * to the box.  The purpose of the event is then to notify listeners to prepare
 * themselves for a different number of molecules.
 */
public interface IBoxMoleculeCountEvent extends IBoxEvent {

    /**
     * @return the species whose number of molecules has changed
     */
    public ISpecies getSpecies();

    /**
     * @return the new number of molecules
     */
    public int getCount();
}
