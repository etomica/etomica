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

}