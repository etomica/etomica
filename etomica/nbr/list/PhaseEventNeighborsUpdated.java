package etomica.nbr.list;

import etomica.phase.Phase;
import etomica.phase.PhaseEvent;

/**
 * Phase event that informs listeners that the neighbor lists in the given
 * Phase have been udpated.
 *
 * @author Andrew Schultz
 */
public class PhaseEventNeighborsUpdated extends PhaseEvent {

    public PhaseEventNeighborsUpdated(Phase phase) {
        super(phase);
    }

    private static final long serialVersionUID = 1L;
}
