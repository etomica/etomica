package etomica.potential;

import etomica.Phase;
import etomica.Space;
import etomica.nbr.NeighborCriterion;

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

    /**
     * Informs the potential of the phase on which it acts. Typically this
     * requires at least that it update the nearestImageTransformer of its
     * coordinatePair (if it uses one), e.g.:
     * cPair.setNearestImageTransformer(phase.boundary());
     */
    public abstract void setPhase(Phase phase);

    public void setCriterion(NeighborCriterion criterion) {
        this.criterion = criterion;
    }

    public NeighborCriterion getCriterion() {
        return criterion;
    }

    protected NeighborCriterion criterion = NeighborCriterion.ALL;
}//end of Potential2

