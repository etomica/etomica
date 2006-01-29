/*
 * Created on Mar 2, 2005
 *
 */
package etomica.nbr;

import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.phase.Phase;

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
public abstract class CriterionAdapter extends NeighborCriterion {

    /**
     * Constructs criterion that wraps the given criterion.
     * 
     * @param criterion
     */
    public CriterionAdapter(NeighborCriterion criterion) {
        subCriterion = criterion;
    }

    /**
     * Implementation of this method should introduce new criterion and return
     * true if pair meets this criterion and that of the wrapped
     * NeighborCriterion.
     */
    public abstract boolean accept(AtomSet pair);

    /**
     * Indicates whether the neighbor list for the given atom should be updated,
     * according to the wrapped criterion.
     * 
     * @return true if the atom's list should be updated.
     */
    public boolean needUpdate(Atom atom) {
        return subCriterion.needUpdate(atom);
    }

    /**
     * Specifies to the wrapped criterion the phase where the criterion is being
     * applied. Sometimes needed if the criterion depends on features of the
     * phase, such as the volume.
     */
    public void setPhase(Phase phase) {
        subCriterion.setPhase(phase);
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
    public void reset(Atom atom) {
        subCriterion.reset(atom);
    }

    /**
     * Returns the nominal distance within which two atoms are considered
     * neighbors, according to the wrapped criterion.
     */
    public boolean isRangeDependent() {
        return subCriterion.isRangeDependent();
    }

    protected final NeighborCriterion subCriterion;
}