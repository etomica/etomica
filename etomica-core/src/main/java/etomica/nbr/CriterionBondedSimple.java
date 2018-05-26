/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtom;

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
    public boolean accept(IAtom atom1, IAtom atom2) {
        int diff = atom1.getIndex() - atom2.getIndex();
        if (isBonded != (diff == 1 || diff == -1) 
                || (atom1.getParentGroup() != atom2.getParentGroup())) {
            return false;
        }
        return subCriterion.accept(atom1, atom2);
    }
    
    private static final long serialVersionUID = 1L;
    private boolean isBonded;
}
