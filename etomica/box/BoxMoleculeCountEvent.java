/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeCountEvent;
import etomica.api.ISpecies;

public class BoxMoleculeCountEvent extends BoxEvent implements IBoxMoleculeCountEvent {

    public BoxMoleculeCountEvent(IBox box, ISpecies _species, int _count) {
        super(box);
        this.species = _species;
        this.count = _count;
    }
    
    public ISpecies getSpecies() {
        return species;
    }
    
    public int getCount() {
        return count;
    }
    
    protected int count = -1;
    protected ISpecies species = null;
}
