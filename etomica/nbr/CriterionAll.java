package etomica.nbr;

import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.phase.Phase;

/**
 * Specifies that all atoms pairs are to be considered neighbors.  Should
 * not be used for species in which atoms are being added/removed by integrator.
 */
public class CriterionAll implements NeighborCriterion, java.io.Serializable {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     * This is appropriate if atoms are never added to or removed from phase,
     * because all atoms are always on neighbor list.
     */
    public boolean needUpdate(IAtom atom) {
        return false;
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
    public void reset(IAtom atom) {
    }

    /**
     * Always returns true, indicating that all atoms pairs are neighbors.
     */
    public boolean accept(AtomSet pair) {
        return true;
    }
    
}
