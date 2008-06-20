package etomica.space3d;

import etomica.api.IVector;
import etomica.space.IOrientation;

/**
 * Interface for an orientation in 3D.  The orientation is defined by a single
 * vector, so the oriented object must have linear symmetry for this to make
 * sense.  IOrientationFull3D can be used for an oriented object without linear
 * symmetry.
 * 
 * @author Andrew Schultz
 */
public interface IOrientation3D extends IOrientation {
    
    /**
     * Rotate the orientation by the amount dt (radians) about the given axis.
     * If the given axis is considered "up" (out of the plane and toward the
     * viewer), then positive rotation is in the counter-clockwise direction.
     */
    public void rotateBy(double dt, IVector axis);
}
