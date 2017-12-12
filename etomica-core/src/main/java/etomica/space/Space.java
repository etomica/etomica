/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;

/**
 * Superclass for all Space instances. Space classes are factories for vectors, tensors, etc.
 */
public abstract class Space implements java.io.Serializable {

    protected Space() {
    }
    
    /**
     * @return a space instance of the given dimension, which must be 1, 2, or 3.
     * @param D the dimension of the requested space instance
     * @throws IllegalArgumentException if D is not 1, 2, or 3.
     */
    public static Space getInstance(int D) {
        switch(D) {
        case 1: 
            return Space1D.getInstance();
        case 2:
            return Space2D.getInstance();
        case 3:
            return Space3D.getInstance();
        default:
            throw new IllegalArgumentException("Illegal dimension for space");
        }
    }

    /**
     * @return the dimension of this space.
     */
    public abstract int D();

    public int getD() {
        return this.D();
    }

    /**
     * @return the given value raised to the power 1/D, where D is the dimension of the space.
     */
    public abstract double rootD(double a);

    /**
     * @return the given value raised to the Dth power, where D is the dimension of the space.
     */
    public abstract int powerD(int a);

    /**
     * @return the given value raised to the Dth power, where D is the dimension of the space.
     */
    public abstract double powerD(double a);

    /**
     * @return a new Vector appropriate to the space.
     */
    public abstract Vector makeVector();

    /**
     * @return a new Orientation appropriate to the space.
     */
    public abstract IOrientation makeOrientation();

    /**
     * @return a new Tensor appropriate to the space.
     */
    public abstract Tensor makeTensor();

    /**
     * @return a new RotationTensor appropriate to the space.
     */
    public abstract RotationTensor makeRotationTensor();

    /**
     * @return an array of dimension D (as defined for this Space instance),
     * with each element equal to the given value.
     */
    public abstract int[] makeArrayD(int i);

    /**
     * @return an array of dimension D (as defined for this Space instance),
     * with each element equal to the given value.
     */
    public abstract double[] makeArrayD(double d);

    /**
     * The "volume" of the "sphere" defined in the D-dimensional space.
     * In 1-D, this is twice the radius; in 2-D the area of the circle;
     * in 3-D the volume of the sphere.
     *
     * @param r the radius
     * @return the volume
     */
    public abstract double sphereVolume(double r);

    /**
     * The surface "area" of the "sphere" defined in the D-dimensional space.
     * In 1-D this is zero; in 2-D the circumference of the circle; in 3-D
     * the surface area of the sphere.
     *
     * @param r the sphere radius
     * @return the area
     */
    public abstract double sphereArea(double r);
    

    /**
     * Returns a Vector from the space of the given dimension.
     * 
     * @throws IllegalArgumentException if D is not 1, 2, or 3.
     */
    public static Vector makeVector(int D) {
        switch(D) {
            case 1:  return new Vector1D();
            case 2:  return new Vector2D();
            case 3:  return new Vector3D();
            default: throw new IllegalArgumentException("Space.makeVector: Requested dimension not implemented");
        }
    }

    /**
     * Returns a Vector initialized to the given set of values in the array.
     * Spatial dimension of the Vector is determined by the length of a.
     * 
     * @throws IllegalArgumentException if a.length is not 1, 2, or 3.
     */
    public Vector makeVector(double[] a) {
        switch(a.length) {
            case 1:  return new etomica.space1d.Vector1D(a);
            case 2:  return new etomica.space2d.Vector2D(a);
            case 3:  return new etomica.space3d.Vector3D(a);
            default: throw new IllegalArgumentException("Space.makeVector: Requested dimension not implemented");
        }
    }
    
    /**
     * Returns a Vector initialized to the given set of values in the array (cast to double).
     * Spatial dimension of the Vector is determined by the length of a.
     */
    public Vector makeVector(int[] k) {
        double[] a = new double[k.length];
        for(int i=0; i<k.length; i++) {a[i] = k[i];}
        return makeVector(a);
    }

    /**
     * Makes an array of Vectors.
     * @param n number of vectors in the returned array
     * @return an array of n new vectors made by this Space instance
     */
    public Vector[] makeVectorArray(int n) {
        Vector[] vectors = new Vector[n];
        for(int i=0; i<n; i++) vectors[i] = makeVector();
        return vectors;
    }
   
}//end of Space    
