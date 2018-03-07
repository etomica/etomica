/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtomList;

/**
 * Pair criterion that judges whether two atoms are or are not in
 * the same molecule. Configurable to accept intra- or inter-molecular
 * pairs. 
 */
public class CriterionMolecular extends CriterionAdapter {

    /**
     * Constructs criterion in default state such that intramolecular pairs are rejected,
     * and intermolecular pairs are accepted.
     */
    public CriterionMolecular(NeighborCriterion criterion) {
        super(criterion);
    }
   
    /**
     * Configures to accept intra- (if argument is true) or inter-
     * (if argument is false) molecular pairs.  Default is false.
     */
    public void setIntraMolecular(boolean b) {
        isIntraMolecular = b;
    }

    /**
     * Flag indicating whether to accept intra- (if argument is true) or inter-
     * (if argument is false) molecular pairs.
     */
    public boolean isIntraMolecular() {
        return isIntraMolecular;
    }
    
    /**
     * Returns false if pair is/isn't in same molecule (depending on setting
     * of intraMolecular); if matches this criterion, return value will be
     * that given by any subCriterion.
     */
    public boolean accept(IAtomList pair) {
        if (isIntraMolecular != (pair.get(0).getParentGroup() == pair.get(1).getParentGroup())) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private static final long serialVersionUID = 1L;
    private boolean isIntraMolecular;
}
