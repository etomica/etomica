package etomica.space;

import etomica.space1d.Space1D;
import etomica.space1d.Vector1D;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;

public abstract class Space implements java.io.Serializable {

    protected Space() {
    }
    
    /**
     * Returns a space instance of the given dimension, which must be 1, 2, or 3.
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
    
    public abstract int D();
    
    public abstract double rootD(double a);
    
    /**
     * Returns the given value raised to the Dth power, where D is the dimension of the space.
     */
    public abstract int powerD(int a);
    /**
     * Returns the given value raised to the Dth power, where D is the dimension of the space.
     */
    public abstract double powerD(double a);
    
    public abstract Vector makeVector();
    
    public abstract Orientation makeOrientation();
    
    public abstract Tensor makeTensor();
    
    public abstract Tensor makeRotationTensor();
    
    /**
     * Returns an array of dimension D, with each element equal to the given value.
     */
    public abstract int[] makeArrayD(int i);
    /**
     * Returns an array of dimension D, with each element equal to the given value.
     */
    public abstract double[] makeArrayD(double d);

    public abstract double sphereVolume(double r);

    public abstract double sphereArea(double r);
    
    /**
     * Returns the square distance between the two vectors, using the given boundary condition.
     */
    public static double r2(Vector u1, Vector u2, Boundary b) { //square distance between two vectors, subject to boundary b
        if(u1.D() != u2.D()) throw new IllegalArgumentException("Space.r2:  Dimension of vectors not equal to each other");
        switch(u1.D()) {
            case 1: return Space1D.r2((Vector1D)u1, (Vector1D)u2, b, new Vector1D());
            case 2: return Space2D.r2((Vector2D)u1, (Vector2D)u2, b, new Vector2D());
            case 3: return Space3D.r2((Vector3D)u1, (Vector3D)u2, b, new Vector3D());
            default: throw new UnsupportedOperationException("Unsupported Space dimension");
        }
    }
    /**
     * Returns a Vector from the space of the given dimension.
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
     */
    public static Vector makeVector(double[] a) {
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
    public static Vector makeVector(int[] k) {
        double[] a = new double[k.length];
        for(int i=0; i<k.length; i++) {a[i] = k[i];}
        return makeVector(a);
    }
    
    /**
     * Instance methods that makes and returns an array of vectors having the
     * given number of elements.
     * @param n number of vectors in the returned array
     * @return an array of n new vectors made by the space instance
     */
    public Vector[] makeVectorArray(int n) {
        Vector[] vectors = new Vector[n];
        for(int i=0; i<n; i++) vectors[i] = makeVector();
        return vectors;
    }
   
}//end of Space    
