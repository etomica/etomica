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

public class OrientationFull3D implements IOrientationFull3D, Serializable {

    /**
     * Default constructor sets orientation to point in the X direction.
     */
    public OrientationFull3D(Space space) {
        direction = space.makeVector();
        direction.setX(0, 1);
        secondaryDirection = space.makeVector();
        secondaryDirection.setX(1, 1);
        v2 = space.makeVector();
        v3 = space.makeVector();
        rotationTensor = (Tensor3D)space.makeTensor();
    }

    public void E(IOrientation o) {
        setDirections(o.getDirection(), ((IOrientationFull3D)o).getSecondaryDirection());
    }
    
    public Vector getDirection() {
        return direction;
    }
    
    public Vector getSecondaryDirection() {
        return secondaryDirection;
    }

    /**
     * Sets this orientation to point in the given direction.
     * @throws Exception if vector has 0 length
     */
    public void setDirection(Vector newDirection) {
        direction.E(newDirection);
        direction.normalize();
    }
    
    public void setDirections(Vector newPrimaryDirection, Vector newSecondaryDirection) {
        direction.E(newPrimaryDirection);
        direction.normalize();
        secondaryDirection.E(newSecondaryDirection);
        secondaryDirection.normalize();
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
        double sindt = Math.sin(dt);
        double cosdt = Math.cos(dt);
        double oneminuscosdt = 1 - cosdt;
        if (Math.abs(dt) < 0.1) {
            // this is a better approximation for small angles (when 1-cos is
            // close to 0)
            oneminuscosdt = sindt*sindt / (1+cosdt);
        }
        // v1 = v1overAxis * axis
        double v1overAxis = axis.dot(direction);
        v3.E(axis);
        v3.XE(direction);

        direction.TE(cosdt);
        direction.PEa1Tv1(oneminuscosdt*v1overAxis, axis);
        direction.PEa1Tv1(sindt, v3);

        // repeat with secondaryDirection
        // v1 = v1overAxis * axis
        v1overAxis = axis.dot(secondaryDirection);
        v3.E(axis);
        v3.XE(secondaryDirection);
        // now temp2 = v3
        secondaryDirection.TE(cosdt);
        secondaryDirection.PEa1Tv1(oneminuscosdt*v1overAxis, axis);
        secondaryDirection.PEa1Tv1(sindt, v3);
    }
    
    /**
     * Applies a random rotation of angle selected uniformly from 0 to tStep.
     * The rotation axis is selected randomly.
     */
    public void randomRotation(IRandom random, double tStep) {
        double tempSq = 0;
        do {
            // first get a random unit vector
            v2.setRandomSphere(random);
            // find the component of the unit vector perpendicular to our direction
            v2.PEa1Tv1(-v2.dot(direction), direction);
            // if the random unit vector was nearly parallel (or anti-parallel)
            // to direction then the calculations will not be particularly
            // precise, so try again
            tempSq = v2.squared();
        } while (tempSq < 0.001);

        v2.TE(1/Math.sqrt(tempSq));
        // we're rotating about the vector perpendicular to both direction and v2
        // (we'll need the axis for secondaryDirection) direction and v2 are both unit vectors
        v3.E(direction);
        v3.XE(v2);
        
        double dt = tStep * random.nextDouble();
        
        // new direction is in the plane of the old direction and temp
        // with components equal to cos(dt) and sin(dt).
        // dt=0 ==> cos(dt)=1 ==> old direction
        double sindt = Math.sin(dt);
        double cosdt = Math.cos(dt);
        direction.TE(cosdt);
        direction.PEa1Tv1(sindt, v2);

        double oneminuscosdt = 1 - cosdt;
        if (Math.abs(dt) < 0.1) {
            oneminuscosdt = sindt*sindt / (1+cosdt);
        }

        v2.E(v3);
        // now v2 is the axis of rotation

        double v1overAxis = v2.dot(secondaryDirection);
        v3.E(v2);
        v3.XE(secondaryDirection);
        // now temp2 = v3
        secondaryDirection.TE(cosdt);
        secondaryDirection.PEa1Tv1(oneminuscosdt*v1overAxis, v2);
        secondaryDirection.PEa1Tv1(sindt, v3);
    }

    /**
     * Normalizes both direction and secondaryDirection
     */
    public void relax() {
        direction.TE(1/Math.sqrt(direction.squared()));
        secondaryDirection.TE(1/Math.sqrt(secondaryDirection.squared()));
    }

    private static final long serialVersionUID = 1L;
    protected final Vector direction, secondaryDirection;
    protected final Vector v2, v3;
    protected final Tensor3D rotationTensor;
}
