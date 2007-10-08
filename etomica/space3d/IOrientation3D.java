package etomica.space3d;

import etomica.space.IVector;
import etomica.space.IOrientation;

/**
 * Interface for a class that specifies an orientation in space.
 */
public interface IOrientation3D extends IOrientation {
    
    /**
     * Rotate the orientation by the amount dt (radians) about the given axis.
     * If the given axis is considered "up" (out of the plane and toward the
     * viewer), then positive rotation is in the counter-clockwise direction.
     */
    public void rotateBy(double dt, IVector axis);
}
