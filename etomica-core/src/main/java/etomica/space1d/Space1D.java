/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space1d;

import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space.RotationTensor;
import etomica.space.Space;

/**
 * Factory and methods appropriate to a 1-dimensional space.  This is
 * a singleton class that can be accessed only via the static getInstance method.
 */
public final class Space1D extends Space {
    
    private static final Space1D INSTANCE = new Space1D();
    private static final long serialVersionUID = 1L;

    /**
     * Private constructor for singleton.
     */
    private Space1D() {
        super();
    }

    /**
     * @return the singleton instance of this class.
     */
    public static Space1D getInstance() {
        return INSTANCE;
    }

    /**
     * Computes the square of the magnitude of the difference of two vectors, subject
     * to a nearest image transformation.  This method constructs a new vector that
     * is used as the work-vector input to the other r2 method.
     */
    public static final double r2(Vector1D u1, Vector1D u2, Boundary b) {
        return r2(u1, u2, b, new Vector1D());
    }
    
    /**
     * Computes the square of the magnitude of the difference of two vectors, subject
     * to a nearest image transformation.
     * @param u1 one of the vectors and is unchanged by the calculation
     * @param u2 the other vector and is unchanged by the calculation
     * @param b a nearest image transformation
     * @param work a work vector used for the calculation
     * @return the nearest-image squared distance
     */
    public static final double r2(Vector1D u1, Vector1D u2, Boundary b, Vector1D work) {
        work.Ev1Mv2(u1, u2);
        b.nearestImage(work);
        return work.squared();
    }

    public int D() {
        return 1;
    }

    public int powerD(int n) {
        return n;
    }

    public double powerD(double a) {
        return a;
    }

    public double rootD(double a) {return a;}

    public double sphereVolume(double r) {
        return 2.0 * r;
    } //volume of a sphere of radius r

    public double sphereArea(double r) {
        return 2.0;
    } //surface area of sphere of radius r (used for differential shell volume)

    public Vector makeVector() {
        return new Vector1D();
    }

    public etomica.space.IOrientation makeOrientation() {
        return new Orientation1D();
    }

    public etomica.space.Tensor makeTensor() {
        return new Tensor1D();
    }

    public RotationTensor makeRotationTensor() {
        return new RotationTensor1D();
    }
    
    public int[] makeArrayD(int i) {
        return new int[] { i };
    }
    
    public double[] makeArrayD(double d) {
        return new double[] { d };
    }

    /**
     * Required to guarantee singleton when deserializing.
     *
     * @return the singleton INSTANCE
     */
    private Object readResolve() {
        return INSTANCE;
    }

}
