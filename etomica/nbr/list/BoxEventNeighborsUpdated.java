package etomica.nbr.list;

import etomica.box.Box;
import etomica.box.BoxEvent;

/**
 * Box event that informs listeners that the neighbor lists in the given
 * Box have been udpated.
 *
 * @author Andrew Schultz
 */
public class BoxEventNeighborsUpdated extends BoxEvent {

    public BoxEventNeighborsUpdated(Box box) {
        super(box);
    }

    private static final long serialVersionUID = 1L;
}
