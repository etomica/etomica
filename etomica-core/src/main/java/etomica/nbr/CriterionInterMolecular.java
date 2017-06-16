/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtomList;

/**
 * Pair criterion that judges whether two atoms are or are not in the same 
 * molecule.  Intermolecular pairs are always accepted.  An optional 
 * intramolecular may be used to filter out intramolecular pairs.
 */
public class CriterionInterMolecular extends CriterionAdapter {

    /**
     * Constructs criterion in default state such that intramolecular pairs are rejected,
     * and intermolecular pairs are accepted.
     */
    public CriterionInterMolecular(NeighborCriterion criterion) {
        super(criterion);
    }
   
    /**
     * Configures to use the given criterion for intramolecular pairs.
     */
    public void setIntraMolecularCriterion(NeighborCriterion criterion) {
        intraCriterion = criterion;
    }

    /**
     * Returns the intramolecular criterion, or null if none is in use.
     */
    public NeighborCriterion getIntraMolecularCriterion() {
        return intraCriterion;
    }
    
    public boolean accept(IAtomList pair) {
        // Only ask the intracriterion if it exists and the pair is intramolecular. 
        if ((pair.getAtom(0).getParentGroup() == pair.getAtom(1).getParentGroup()) && (intraCriterion == null ||
                !intraCriterion.accept(pair))) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private static final long serialVersionUID = 1L;
    private NeighborCriterion intraCriterion;
}
