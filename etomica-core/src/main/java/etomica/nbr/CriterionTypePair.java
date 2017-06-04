/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.nbr;

import etomica.atom.AtomType;
import etomica.atom.IAtomList;

/**
 * Filters atoms pairs to match a given pair of AtomTypes.
 * 
 * @author Andrew Schultz
 */
public class CriterionTypePair extends CriterionAdapter {

    private static final long serialVersionUID = 1L;
    private final AtomType type0;
    private final AtomType type1;

    public CriterionTypePair(NeighborCriterion criterion,
                             AtomType type0, AtomType type1) {
        super(criterion);
        this.type0 = type0;
        this.type1 = type1;
    }

    /**
     * Returns true if the AtomTypes for the pair of atoms match the AtomTypes
     * given at construction (without regard to the order of the pair), and if
     * the wrapped criterion also accepts the pair.
     */
    public boolean accept(IAtomList pair) {
        AtomType atom0Type = pair.getAtom(0).getType();
        AtomType atom1Type = pair.getAtom(1).getType();
        if ( (atom0Type == type0 && atom1Type == type1) ||
             (atom0Type == type1 && atom1Type == type0) ) {
            return subCriterion.accept(pair);
        }
        return false;
    }

    /**
     * Returns the AtomTypes accepted by this NeighborCriterion
     */
    public AtomType[] getTypes() {
        return new AtomType[]{type0, type1};
    }
}
