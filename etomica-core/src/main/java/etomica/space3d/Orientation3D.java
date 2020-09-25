/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space3d;

import java.io.Serializable;

import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.util.Debug;

public class Orientation3D implements IOrientation3D, Serializable {

    private static final long serialVersionUID = 1L;
    protected final Vector direction;
    protected final Vector temp, temp2;

    /**
     * Default constructor sets orientation to point in the X direction.
     */
    public Orientation3D(Space space) {
        this(space.makeVector());
    }

    public Orientation3D(Vector direction) {
        this.direction = direction;
        this.direction.setX(0, 1);
        temp = Space.makeVector(3);
        temp2 = Space.makeVector(3);
    }

    public void E(IOrientation o) {
        setDirection(o.getDirection());
    }

    public Vector getDirection() {
        return direction;
    }

    /**
     * Sets this orientation to point in the given direction.
     * @throws Exception if vector has 0 length
     */
    public void setDirection(Vector newDirection) {
        direction.E(newDirection);
        direction.normalize();
    }

    /**
     * Rotates orientation by the given value about the given axis.  The axis
     * must have unit length, but need not be perpendicular to the current
     * orientation direction.
     */
    public void rotateBy(double dt, Vector axis) {
        // consider a circle on the surface of the unit sphere.  The given axis
        // passes through the center of the circle.  The circle passes through
        // the current direction vector and the vector v4 defined below.  We
        // rotate the direction by the given angle (dt) around the circle.

        // v1 is the projection of direction onto axis
        // v2 is the component of direction perpendicular to axis
        // v3 has the same magnitude as v2 and is perpendicular to both
        //    direction and axis
        // v4 is a unit vector whose components are v1 and v3
        if (Debug.ON && Math.abs(axis.squared() - 1) > 1E-10) {
            throw new IllegalArgumentException("I need a unit vector for the axis");
        }
        // v1 = v1overAxis * axis
        double v1overAxis = axis.dot(direction);

        temp.Ea1Tv1(-v1overAxis, axis);
        temp.PE(direction);
        // now temp = v2
        temp2.E(axis);
        temp2.XE(direction);
        // now temp2 = v3
        direction.Ea1Tv1(Math.cos(dt), temp);
        direction.PEa1Tv1(Math.sin(dt), temp2);
        direction.PEa1Tv1(v1overAxis, axis);
    }

    /**
     * Applies a random rotation of angle selected uniformly from 0 to tStep.
     * The rotation axis is selected randomly, but is orthogonal to the
     * original direction.
     */
    public void randomRotation(IRandom random, double tStep) {
        double tempSq = 0;
        do {
            // first get a random unit vector
            temp.setRandomSphere(random);
            // find the component of the unit vector perpendicular to our direction
            temp.PEa1Tv1(-temp.dot(direction), direction);
            // if the random unit vector was nearly parallel (or anti-parallel)
            // to direction then the calculations will not be particularly
            // precise, so try again
            tempSq = temp.squared();
        } while (tempSq < 0.001);

        temp.TE(1/Math.sqrt(tempSq));
        double dt = tStep * random.nextDouble();

        // new direction is in the plane of the old direction and temp
        // with components equal to cos(dt) and sin(dt).
        // dt=0 ==> cos(dt)=1 ==> old direction
        direction.TE(Math.cos(dt));
        direction.PEa1Tv1(Math.sin(dt), temp);
    }
}
