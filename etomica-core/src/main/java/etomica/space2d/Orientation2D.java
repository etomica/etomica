/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space2d;

import java.io.Serializable;

import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.util.Constants;

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
     * @throws Exception if vector has 0 length
     */
    public Orientation2D(Vector direction) {
        this.direction = Space.makeVector(2);
        setDirection(direction);
    }
    
    public void E(IOrientation o) {
        setDirection(o.getDirection());
    }
    
    public Vector getDirection() {
        return direction;
    }
    
    public void setDirection(Vector newDirection) {
        direction.E(newDirection);
        direction.normalize();
        angle = Math.atan2(this.direction.getX(1), this.direction.getX(0));
    }
    
    /**
     * Rotates orientation around the given axis by the given value.
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
    protected Vector direction;
}
