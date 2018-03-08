/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.nbr;

import etomica.atom.AtomType;
import etomica.atom.IAtom;

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
     * @param atom1
     * @param atom2
     */
    public boolean accept(IAtom atom1, IAtom atom2) {
        return atom1.getType() == type && subCriterion.accept(atom1, atom2);
    }

    /**
     * Returns the AtomType accepted by this criterion.
     */
    public AtomType getType() {
        return type;
    }
}
