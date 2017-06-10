/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.geometry;

import etomica.space.Vector;
import etomica.exception.MethodNotImplementedException;


/**
 * An arbitrary-dimension sphere.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 28, 2005 by kofke
 */
public class Hypersphere implements Shape, java.io.Serializable {

    /**
     * Creates hypersphere of the given dimension, positioned at the origin.
     */
    public Hypersphere(int D) {
        this(D, 1.0);
    }

    /**
     * Creates hypersphere of the given dimension and radius, positioned at the origin.
     */
    public Hypersphere(int D, double radius) {
        this.D = D;
        setRadius(radius);
    }
    
    public double getVolume() {
        //TODO implement this
        throw new MethodNotImplementedException();
    }

    /* (non-Javadoc)
     * @see etomica.math.geometry.Shape#D()
     */
    public final int D() {
        return D;
    }

    /**
     * Returns true if the given point lies on or in the sphere.
     * 
     * @throws ClassCastException
     *             if given vector is not an instance of a D-dimensional vector
     */
    public boolean contains(Vector r) {
        return r.Mv1Squared(position) <= radius;
    }

    /**
     * @return Returns the radius.
     */
    public double getRadius() {
        return radius;
    }

    /**
     * @param radius
     *            The radius to set.
     */
    public void setRadius(double radius) {
        this.radius = radius;
    }

    /**
     * @return Returns the vector (not a copy) used to represent the position of
     *         the sphere.
     */
    public Vector getPosition() {
        return position;
    }

    /**
     * @param position
     *            The new position, which is copied to an internal vector.
     */
    public void setPosition(Vector position) {
        this.position.E(position);
    }

    protected double radius;
    private Vector position;
    private final int D;

}
