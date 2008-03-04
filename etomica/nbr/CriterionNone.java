package etomica.nbr;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IBox;

public final class CriterionNone implements NeighborCriterion, java.io.Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     * This is appropriate if atoms are never added to or removed from box,
     * because all atoms are always on neighbor list.
     */
    public boolean needUpdate(IAtom atom) {return false;}

    /**
     * Performs no action.
     */
    public void setBox(IBox box) {}

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     */
    public boolean unsafe() {return false;}

    /**
     * Performs no action.
     */
    public void reset(IAtom atom) {}

    /**
     * Always returns false, indicating that no atoms pairs are neighbors.
     */
    public boolean accept(IAtomSet pair) {return false;}
}