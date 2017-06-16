/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.nbr;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;

/**
 * Filters atoms to match a given AtomType.
 * 
 * @author Andrew Schultz
 */
public class CriterionType extends CriterionAdapter {

    private static final long serialVersionUID = 1L;
    private final AtomType type;

    public CriterionType(NeighborCriterion criterion,
                         AtomType type) {
        super(criterion);
        this.type = type;
    }
    
    /**
     * Returns true if the AtomType of the atom matches the AtomType given at
     * construction and if the wrapped criterion accept also returns true.
     */
    public boolean accept(IAtomList atom) {
        if (atom.getAtom(0).getType() == type) {
            return subCriterion.accept(atom);
        }
        return false;
    }

    /**
     * Returns the AtomType accepted by this criterion.
     */
    public AtomType getType() {
        return type;
    }
}
