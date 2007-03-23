package etomica.space;

import etomica.util.IRandom;

/**
 * Interface for a class that specifies an orientation in space.
 */
//TODO consider if this should be an interface
public abstract class Orientation implements java.io.Serializable {
    
    /**
     * Copies the given orientation to this one.
     */
    public abstract void E(Orientation o);

    /**
     * Set of angles describing the orientation, according to a
     * convention defined by the implementation.
     */
    public abstract double[] angle();

    /**
     * Rotate all orientation angles by the corresponding amounts in the
     * given array.
     */
    public abstract void rotateBy(double[] t);

    /**
     * Rotate the angle indexed i by the amount dt.
     */
    public abstract void rotateBy(int i, double dt);

    /**
     * Perform a rotation by a random amount in the solid angle theta on the 
     * present orientation.
     */
    public abstract void randomRotation(IRandom random, double theta);

    /**
     * Changes the components of all vectors in the given array 
     * to the body frame representation for the current orientation.
     */
    public abstract void convertToBodyFrame(IVector[] v);

    /**
     * Changes the components all vectors in the given array from a body frame of
     * the current orientation to the space frame.
     */
    public abstract void convertToSpaceFrame(IVector[] v);
}
