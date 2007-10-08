package etomica.space;

import etomica.util.IRandom;

/**
 * Interface for a class that specifies an orientation in space.
 */
public interface IOrientation {
    
    /**
     * Copies the given orientation to this one.
     */
    public void E(IOrientation o);

    /**
     * Returns a unit vector pointing in the orientation's direction.  This
     * vector should not be modified.
     */
    public IVector getDirection();
    
    /**
     * Sets the orientation's direction to be the given direction.
     */
    public void setDirection(IVector newDirection);
    
    /**
     * Perform a rotation by a random amount in the solid angle theta on the 
     * present orientation.
     */
    public void randomRotation(IRandom random, double theta);
}
