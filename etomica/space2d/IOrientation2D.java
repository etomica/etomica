package etomica.space2d;

import etomica.space.IOrientation;

/**
 * Interface for a class that specifies an orientation in space.
 */
public interface IOrientation2D extends IOrientation {
    
    /**
     * Rotate the orientation by the amount dt (radians) in the
     * counter-clockwise direction.
     */
    public void rotateBy(double dt);
}
