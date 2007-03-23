package etomica.space2d;

import etomica.space.IVector;
import etomica.space.Orientation;
import etomica.units.Angle;
import etomica.units.Dimension;
import etomica.util.Constants;
import etomica.util.IRandom;

/**
 * Orientation in a 2-dimensional space. Orientation can be expressed in terms
 * of a single angle that describes rotation of the body frame in a
 * counter-clockwise direction with respect to the fixed space frame.  Periodicity is
 * applied to restrict the angle to values between 0 and 2 PI.
 */
public class Orientation2D extends Orientation {

    /**
     * Default constructor sets orientation angle to zero.
     */
    public Orientation2D() {
        this(0.0);
    }

    /**
     * Constructs with orientation as specified by the given angle theta.
     */
    public Orientation2D(double theta) {
        setTheta(theta);
    }
    
    /**
     * Sets this orientation to that specified by the given orientation.
     * 
     * @throws ClassCastException
     *             if the given orientation is not of type Orientation2D
     */
    public void E(Orientation o) {
        setTheta(((Orientation2D) o).angle[0]);
    }

    /**
     * Returns theta.
     */
    public double[] angle() {
        return angle;
    }

    /**
     * Rotates orientation by the given value.
     * 
     * @throws IllegalArgumentException if index is not zer0.
     */
    public void rotateBy(int index, double dt) {
        if (index == 0) {
            rotateBy(dt);
        } else {
            throw new IllegalArgumentException("Invalid angle index: "+index);
        }
    }

    /**
     * Rotates orientation the value of the first element of the array.
     */
    public void rotateBy(double[] dt) {
        rotateBy(dt[0]);
    }

    /**
     * Rotates angle by the given angle.  Rotation is counter-clockwise
     * if angle is positive, clockwise if it is negative.  
     */
    public void rotateBy(double dt) {
        setTheta(angle[0]+dt);
    }
    
    /**
     * Sets orientation to that of the given angle.  
     */
    public void setTheta(double theta) {
        while (theta > Constants.TWO_PI) {
            theta -= Constants.TWO_PI;
        }
        while (theta < 0.0) {
            theta += Constants.TWO_PI;
        }
        angle[0] = theta;
    }
    
    public double getTheta() {
        return angle[0];
    }
    
    /**
     * Returns Angle.DIMENSION, indicating that theta is an angle.
     */
    public Dimension getThetaDimension() {
        return Angle.DIMENSION;
    }

    /**
     * Applies a random rotation of angle selected uniformly within plus or minus tStep 
     * from the current value.
     */
    public void randomRotation(IRandom random, double tStep) {
        rotateBy((2. * random.nextDouble() - 1.0) * tStep);
    }

    /**
     * Takes vectors defined in the space-frame representation and converts
     * each to the equivalent vector defined in the body-frame representation.
     */
    public void convertToBodyFrame(IVector[] v) {
        double axx = Math.cos(angle[0]);
        double axy = -Math.sin(angle[0]);
        for(int i=0; i<v.length; i++) {
            convert((Vector2D)v[i], axx, axy);
        }
    }

    /**
     * Takes vectors defined in the body-frame representation and converts
     * each to the equivalent vector defined in the space-frame representation.
     */
    public void convertToSpaceFrame(IVector[] v) {
        double axx = Math.cos(angle[0]);
        double axy = Math.sin(angle[0]);
        for(int i=0; i<v.length; i++) {
            convert((Vector2D)v[i], axx, -axy);
        }
    }
    
    private void convert(Vector2D v, double axx, double axy) {
        double x = axx * v.x + axy * v.y;
        v.y = -axy * v.x + axx * v.y;
        v.x = x;
        
    }

    private static final long serialVersionUID = 1L;
    private final double[] angle = new double[1];
}
