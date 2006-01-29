/**
 * 
 */
package etomica.nbr;

import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.phase.Phase;

public final class CriterionNone extends NeighborCriterion {
    /**
     * Always returns false, indicating that neighbor list never needs updating.
     * This is appropriate if atoms are never added to or removed from phase,
     * because all atoms are always on neighbor list.
     */
    public boolean needUpdate(Atom atom) {return false;}

    public boolean isRangeDependent() {return false;}

    /**
     * Performs no action.
     */
    public void setPhase(Phase phase) {}

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     */
    public boolean unsafe() {return false;}

    /**
     * Performs no action.
     */
    public void reset(Atom atom) {}

    /**
     * Always returns false, indicating that no atoms pairs are neighbors.
     */
    public boolean accept(AtomSet pair) {return false;}
}