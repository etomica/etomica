package etomica;

import etomica.space.Boundary;
import etomica.space.Coordinate;
import etomica.space.CoordinatePair;
import etomica.space.Orientation;
import etomica.space.Tensor;
import etomica.space.Vector;

/* History of changes
 * 09/01/02 (DAK) added accelerateTo method to Coordinate
 * 07/10/03 (DAK) added resetV method to CoordinatePair
 * 08/27/03 (DAK) added isZero method to Vector
 * 12/09/03 (DAK) added setRandomInSphere method in Vector
 * 01/22/04 (DAK) added (then removed) CoordinateGroup interface
 */

public abstract class Space implements java.io.Serializable {
    
    public final int D;
    private final double rD; // reciprocal of D
    
    protected Space(int d) {
        if(d < 1 || d > 3) throw new IllegalArgumentException("Illegal dimension for space");
        D = d;
        rD = 1.0/D;
    }
    
    /**
     * Returns a space instance of the given dimension, which must be 1, 2, or 3.
     */
    public static Space makeSpace(int D) {
        switch(D) {
        case 1: 
            return Space1D.INSTANCE;
        case 2:
            return Space2D.INSTANCE;
        case 3:
            return Space3D.INSTANCE;
        default:
            throw new IllegalArgumentException("Illegal dimension for space");
        }
    }
    
    public int D() {return D;}
    /**
     * Returns the given value raised to the Dth power, where D is the dimension of the space.
     */
    public abstract int powerD(int a);
    /**
     * Returns the given value raised to the Dth power, where D is the dimension of the space.
     */
    public abstract double powerD(double a);
    
    /**
     * Returns the Dth root of the given value, a^(1/D), where D is the dimension of the space.
     */
    public double rootD(double a) {return Math.pow(a, rD);}
    
    public abstract Vector origin();
    
    public abstract Vector makeVector();      //Space.Vector
    public abstract Orientation makeOrientation();
    public abstract Tensor makeTensor();
    public abstract Tensor makeRotationTensor();
    public abstract Coordinate makeCoordinate(Atom a);
    public abstract CoordinatePair makeCoordinatePair();
    public abstract Boundary makeBoundary();  //makes boundary of default type
    public abstract Boundary makeBoundary(Boundary.Type type);
    public abstract Boundary.Type[] boundaryTypes();
    public boolean requiresSpecialBoundary() {return false;}
    /**
     * Returns an array of dimension D, with each element equal to the given value.
     */
    public abstract int[] makeArrayD(int i);
    /**
     * Returns an array of dimension D, with each element equal to the given value.
     */
    public abstract double[] makeArrayD(double d);

    public double sphereVolume(double r) {
        switch(D) {
            case 1: return 2*r;
            case 2: return Math.PI*r*r;
            default:
            case 3: return (Math.PI*4.0*r*r*r/3.0);
        }
    }
    public double sphereArea(double r)  {
        switch(D) {
            case 1: return 2;
            case 2: return 2*Math.PI*r;
            default:
            case 3: return 4*Math.PI*r*r;
        }
    }
    
    /**
     * Returns the square distance between the two vectors, using the given boundary condition.
     */
    public static double r2(Vector u1, Vector u2, Boundary b) { //square distance between two vectors, subject to boundary b
        if(u1.D() != u2.D()) throw new IllegalArgumentException("Space.r2:  Dimension of vectors not equal to each other");
        switch(u1.D()) {
            case 1: return Space1D.r2((etomica.space1d.Vector)u1, (etomica.space1d.Vector)u2, (Boundary)b);
            case 2: return Space2D.r2((etomica.space2d.Vector)u1, (etomica.space2d.Vector)u2, (Boundary)b);
            case 3: return Space3D.r2((etomica.space3d.Vector)u1, (etomica.space3d.Vector)u2, (Boundary)b);
            default: throw new IllegalArgumentException("Space.r2: Unknown vector dimension");
        }
    }
    /**
     * Returns a Vector from the space of the given dimension.
     */
    public static Vector makeVector(int D) {
        switch(D) {
            case 1:  return new etomica.space1d.Vector();
            case 2:  return new etomica.space2d.Vector();
            case 3:  return new etomica.space3d.Vector();
            default: throw new IllegalArgumentException("Space.makeVector: Requested dimension not implemented");
        }
    }
    /**
     * Returns a Vector initialized to the given set of values in the array.
     * Spatial dimension of the Vector is determined by the length of a.
     */
    public static Vector makeVector(double[] a) {
        switch(a.length) {
            case 1:  return new etomica.space1d.Vector(a);
            case 2:  return new etomica.space2d.Vector(a);
            case 3:  return new etomica.space3d.Vector(a);
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
