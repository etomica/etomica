package etomica.potential;

import etomica.nbr.CriterionAll;
import etomica.nbr.NeighborCriterion;
import etomica.space.Space;

/**
 * Potential acting on or within an atom, or between a pair of atoms or atom
 * groups. Contains other Potential instances that describe the specific
 * interactions between the atoms of the group(s).
 * 
 * @author David Kofke
 */

public abstract class Potential2 extends Potential {

    /**
     * Constructs potential with given space a a default CoordinatePair for
     * spherical, non-kinetic molecules.
     */
    public Potential2(Space space) {
        super(2, space);
    }

    public void setCriterion(NeighborCriterion criterion) {
        this.criterion = criterion;
    }

    public NeighborCriterion getCriterion() {
        return criterion;
    }

    protected NeighborCriterion criterion = new CriterionAll();
}
