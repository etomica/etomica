package etomica.box;

import etomica.api.IBox;
import etomica.api.IVector;

/**
 * Interface for an object which returns random positions within a box's
 * boundary.
 *
 * @author Andrew Schultz
 */
public interface RandomPositionSource {

    /**
     * Notifies the RandomPositionSource of the box from which random positions
     * should be taken from.
     */
    public void setBox(IBox box);
    
    /**
     * Returns a random position with the previously set box's boundary.
     */
    public IVector randomPosition();
}
