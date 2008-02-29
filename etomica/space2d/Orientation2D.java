package etomica.space2d;

import java.io.Serializable;

import etomica.api.IVector;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.util.Constants;
import etomica.util.IRandom;

/**
 * Orientation in a 2-dimensional space. Orientation can be expressed in terms
 * of a single angle that describes rotation of the body frame in a
 * counter-clockwise direction with respect to the fixed space frame.  Periodicity is
 * applied to restrict the angle to values between 0 and 2 PI.
 */
public class Orientation2D implements IOrientation2D, Serializable {

    /**
     * Default constructor sets orientation angle to zero.
     */
    public Orientation2D() {
        direction = Space.makeVector(2);
        direction.setX(0, 1);
        angle = 0;
    }

    /**
     * Constructs with orientation as specified by the given angle theta.
     * @throws an exception if vector has 0 length
     */
    public Orientation2D(IVector direction) {
        this.direction = Space.makeVector(2);
        setDirection(direction);
    }
    
    public void E(IOrientation o) {
        setDirection(o.getDirection());
    }
    
    public IVector getDirection() {
        return direction;
    }
    
    /**
     * Sets this orientation to point in the given direction.
     * @throws an exception if vector has 0 length
     */
    public void setDirection(IVector newDirection) {
        direction.E(newDirection);
        direction.normalize();
        angle = Math.atan2(this.direction.x(1), this.direction.x(0));
    }
    
    /**
     * Rotates orientation around the given axis by the given value.
     * @throws IllegalArgumentException if index is not zer0.
     */
    public void rotateBy(double dt) {
        angle += dt;
        while (angle > Constants.TWO_PI) {
            angle -= Constants.TWO_PI;
        }
        while (angle < 0.0) {
            angle += Constants.TWO_PI;
        }
        direction.setX(0, Math.cos(angle));
        direction.setX(1, Math.cos(angle));
    }

    /**
     * Applies a random rotation of angle selected uniformly within plus or minus tStep 
     * from the current value.
     */
    public void randomRotation(IRandom random, double tStep) {
        rotateBy((2. * random.nextDouble() - 1.0) * tStep);
    }

    private static final long serialVersionUID = 1L;
    protected double angle;
    protected IVector direction;
}
