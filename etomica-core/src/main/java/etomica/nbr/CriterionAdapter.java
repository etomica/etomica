/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on Mar 2, 2005
 *
 */
package etomica.nbr;

import etomica.atom.IAtom;
import etomica.box.Box;

/**
 * Wraps another criterion while adding additional criteria to the acceptance
 * decision. Used to introduce criteria related to the species, molecule, and
 * other identifiers of the pair. Does not normally consider atom positions, but
 * defers this to the wrapped iterator. Consequently the methods needUpdate,
 * unsafe, and others that related to the atom configurations are simply passed
 * on to the wrapped criterion.
 * 
 * @author andrew
 */
public abstract class CriterionAdapter implements NeighborCriterion, java.io.Serializable {

    /**
     * Constructs criterion that wraps the given criterion.
     * 
     * @param criterion
     */
    public CriterionAdapter(NeighborCriterion criterion) {
        subCriterion = criterion;
    }

    /**
     * Returns the criterion wrapped by this adapter
     */
    public final NeighborCriterion getWrappedCriterion() {
        return subCriterion;
    }
    
    /**
     * Implementation of this method should introduce new criterion and return
     * true if pair meets this criterion and that of the wrapped
     * NeighborCriterion.
     * @param atom1
     * @param atom2
     */
    public abstract boolean accept(IAtom atom1, IAtom atom2);

    /**
     * Indicates whether the neighbor list for the given atom should be updated,
     * according to the wrapped criterion.
     * 
     * @return true if the atom's list should be updated.
     */
    public boolean needUpdate(IAtom atom) {
        return subCriterion.needUpdate(atom);
    }

    /**
     * Specifies to the wrapped criterion the box where the criterion is being
     * applied. Sometimes needed if the criterion depends on features of the
     * box, such as the volume.
     */
    public void setBox(Box box) {
        subCriterion.setBox(box);
    }

    /**
     * Indicates whether the atom has changed (e.g. moved) by an amount that
     * might have caused its wrapped-criterion neighbor list to be invalid. If
     * this method returns true, a neighbor list failure may have introduced
     * errors in the calculation.
     */
    public boolean unsafe() {
        return subCriterion.unsafe();
    }

    /**
     * Indicates to wrapped criterion that given atom's neighbor list has just
     * been updated, and that properties (e.g., record of atom's position) used
     * by needUpdate and unsafe() methods should be reset.
     */
    public void reset(IAtom atom) {
        subCriterion.reset(atom);
    }

    protected final NeighborCriterion subCriterion;
}
