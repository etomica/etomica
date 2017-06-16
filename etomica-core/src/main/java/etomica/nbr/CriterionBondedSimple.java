/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtomList;

/**
 * @author andrew
 * Used for bonding potentials. The accept method returns true if the atoms are
 * adjacent in the list of atoms and they have the same parent. 
 */
public class CriterionBondedSimple extends CriterionAdapter {

    public CriterionBondedSimple(NeighborCriterion criterion) {
        super(criterion);
    }
    
    public void setBonded(boolean b) {
        isBonded = b;
    }
    public boolean isBonded() {
        return isBonded;
    }
    
    // always enforce intramolecularity
    public boolean accept(IAtomList pair) {
        int diff = pair.getAtom(0).getIndex() - pair.getAtom(1).getIndex();
        if (isBonded != (diff == 1 || diff == -1) 
                || (pair.getAtom(0).getParentGroup() != pair.getAtom(1).getParentGroup())) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private static final long serialVersionUID = 1L;
    private boolean isBonded;
}
