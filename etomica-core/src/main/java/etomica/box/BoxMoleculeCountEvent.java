/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.ISpecies;

/**
 * Box event that indicates the number of molecules of a particular species has
 * changed.  This event is fired before the molecules have actually been added
 * to the box.  The purpose of the event is then to notify listeners to prepare
 * themselves for a different number of molecules.
 */
public class BoxMoleculeCountEvent extends BoxEvent {

    protected final int count;
    protected final ISpecies species;

    public BoxMoleculeCountEvent(Box box, ISpecies species, int count) {
        super(box);
        this.species = species;
        this.count = count;
    }

    /**
     * @return the species whose number of molecules has changed
     */
    public ISpecies getSpecies() {
        return species;
    }

    /**
     * @return the new number of molecules
     */
    public int getCount() {
        return count;
    }
}
