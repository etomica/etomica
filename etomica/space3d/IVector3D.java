package etomica.space3d;

import etomica.space.IVector;

/**
 * Interface for a 3D vector, which has a XE (cross-product) method
 *
 * @author Andrew Schultz
 */
public interface IVector3D extends IVector {

    /**
     * Sets this vector equal to the cross product of this vector with the
     * given vector.
     */
    public void XE(IVector3D u);

    /**
     * Sets the vector components equal to the given values.
     */
    public void E(double a, double b, double c);
}