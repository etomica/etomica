package etomica.space3d;

import etomica.api.IVector;

public interface IOrientationFull3D extends IOrientation3D {

    /**
     * Returns a unit vector pointing in the orientation's secondary direction.
     * This vector should not be modified.
     */
    public IVector getSecondaryDirection();
    
    /**
     * Sets the orientation's primary and secondary direction to be the given
     * directions.  The two vectors should be orthogonal.
     */
    public void setDirections(IVector newPrimaryDirection, IVector newSecondaryDirection);
}
