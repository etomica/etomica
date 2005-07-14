package etomica.potential;

import etomica.Phase;
import etomica.Potential;
import etomica.Space;
import etomica.nbr.NeighborCriterion;
import etomica.space.CoordinatePair;

/**
 * Potential acting on or within an atom, or between a pair of atoms or atom
 * groups. Contains other Potential instances that describe the specific
 * interactions between the atoms of the group(s).
 * 
 * @author David Kofke
 */

public abstract class Potential2 extends Potential {

    //TODO push coordinate pair into subclasses
    
    /**
     * Constructs potential with given space a a default CoordinatePair
     * for spherical, non-kinetic molecules.
     */
    public Potential2(Space space) {
        this(space, new CoordinatePair(space));
    }

    /**
     * Constructs potential with given space and CoordinatePair.
     */
    public Potential2(Space space, CoordinatePair cPair) {
        super(2, space);
        this.cPair = cPair;
    }

    public void setPhase(Phase phase) {
        cPair.setNearestImageTransformer(phase.boundary());
    }

    public void setCriterion(NeighborCriterion criterion) {
        this.criterion = criterion;
    }

    public NeighborCriterion getCriterion() {
        return criterion;
    }

    protected final CoordinatePair cPair;
    protected NeighborCriterion criterion = NeighborCriterion.ALL;
    
}//end of Potential2

