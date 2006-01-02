package etomica.space;

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
     * Axes defining the orientation of the frame of the oriented body.
     */
    public abstract Vector[] bodyFrame();// body-frame axes in the
                                            // space-fixed frame

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
    public abstract void randomRotation(double theta);

    /**
     * Changes the components of the given vector from a space frame
     * to the body frame representation for the current orientation.
     */
    public abstract void convertToBodyFrame(Vector v);

    /**
     * Changes the components of the given vector from a body frame of
     * the current orientation to the space frame.
     */
    public abstract void convertToSpaceFrame(Vector v);
}
