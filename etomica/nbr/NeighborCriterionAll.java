/*
 * History
 * Created on Nov 27, 2004 by kofke
 */
package etomica.nbr;

import etomica.Atom;
import etomica.AtomSet;
import etomica.Phase;

/**
 * Specifies that all atoms pairs are to be considered neighbors.  Should
 * not be used for species in which atoms are being added/removed by integrator.
 */
public class NeighborCriterionAll extends NeighborCriterion {

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     * This is appropriate if atoms are never added to or removed from phase,
     * because all atoms are always on neighbor list.
     */
    public boolean needUpdate(Atom atom) {
        return false;
    }

    public double getNeighborRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    /**
     * Performs no action.
     */
    public void setPhase(Phase phase) {
    }

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     * This is appropriate if atoms are never added to or removed from phase,
     * because all atoms are always on neighbor list.
     */
    public boolean unsafe() {
        return false;
    }

    /**
     * Performs no action.
     */
    public void reset(Atom atom) {
    }

    /**
     * Always returns true, indicating that all atoms pairs are neighbors.
     */
    public boolean accept(AtomSet a) {
        return true;
    }
    
}
