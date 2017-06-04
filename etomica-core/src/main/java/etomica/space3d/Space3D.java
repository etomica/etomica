/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space3d;

import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space.RotationTensor;
import etomica.space.Space;

/**
 * Factory and methods appropriate to a 3-dimensional space.  This is
 * a singleton class that can be accessed only via the static getInstance method.
 */
public final class Space3D extends Space {

    
    /**
     * Private constructor for singleton.
     */
    private Space3D() {
        super();
    }

    /**
     * @return Returns the instance.
     */
    public static Space3D getInstance() {
        return INSTANCE;
    }

    public final int D() {
        return 3;
    }

    public final int powerD(int n) {
        return n * n * n;
    }

    public final double powerD(double a) {
        return a * a * a;
    }
    
    /**
     * Returns the cube root of the given value, a^(1/D) which is a^(1/3).
     */
    public double rootD(double a) {return Math.pow(a, 1./3.);}
    

    public int[] makeArrayD(int i) {
        return new int[] { i, i, i };
    }

    public double[] makeArrayD(double d) {
        return new double[] { d, d, d };
    }

   public double sphereVolume(double r) {
        return (Math.PI * 4.0 * r * r * r / 3.0);
    }

    public double sphereArea(double r) {
        return (Math.PI * 4 * r * r);
    }

    public Vector makeVector() {
        return new Vector3D();
    }

    public etomica.space.IOrientation makeOrientation() {
        return new OrientationFull3D(this);
    }

    public etomica.space.Tensor makeTensor() {
        return new Tensor3D();
    }

    public RotationTensor makeRotationTensor() {
        return new RotationTensor3D();
    }

    /**
     * Computes the square of the magnitude of the difference of two vectors, subject
     * to a nearest image transformation.  This method constructs a new vector that
     * is used as the work-vector input to the other r2 method.
     */
    public static final double r2(Vector3D u1, Vector3D u2, Boundary b) {
        return r2(u1, u2, b, new Vector3D());
    }

    /**
     * Computes the square of the magnitude of the difference of two vectors, subject
     * to a nearest image transformation.  
     * @param u1 one of the vectors and is unchanged by the calculation
     * @param u2 the other vector and is unchanged by the calculation
     * @param b a nearest image transformation
     * @param work a work vector used for the calculation.
     */
    public static final double r2(Vector3D u1, Vector3D u2, Boundary b,
            Vector3D work) {
        work.Ev1Mv2(u1, u2);
        b.nearestImage(work);
        return work.squared();
    }
    
    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton INSTANCE
     */
    private Object readResolve() {
        return INSTANCE;
    }

    private static final Space3D INSTANCE = new Space3D();
    private static final long serialVersionUID = 1L;

}
