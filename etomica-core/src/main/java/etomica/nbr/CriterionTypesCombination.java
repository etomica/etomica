/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.nbr;

import etomica.atom.IAtomList;
import etomica.atom.IAtomType;

/**
 * Filters AtomSets to match a given set of AtomTypes.  CriterionType and
 * CriterionTypePair should be used for single and pairs of Atoms.
 * 
 * @author Andrew Schultz
 */
public class CriterionTypesCombination extends CriterionAdapter {

    public CriterionTypesCombination(NeighborCriterion criterion, 
            IAtomType[] types) {
        super(criterion);
        this.types = types.clone();
    }
    
    /**
     * Returns true if the AtomTypes for the pair of atoms match the AtomTypes 
     * given at construction (without regard to the order of the AtomSet), and 
     * if the wrapped criterion also accepts the AtomSet.
     */
    public boolean accept(IAtomList atoms) {
        final int nAtoms = atoms.getAtomCount();
        for (int i=0; i<nAtoms; i++) {
            boolean accepted = false;
            for (int j=0; j<types.length; j++) {
                if (atoms.getAtom(i).getType() == types[j]) {
                    accepted = true;
                }
            }
            if (!accepted) {
                return false;
            }
        }
        return subCriterion.accept(atoms);
    }
    
    /**
     * Returns the AtomTypes accepted by this NeighborCriterion
     */
    public IAtomType[] getTypes() {
        return types;
    }
    
    private static final long serialVersionUID = 1L;
    private final IAtomType[] types;
}
