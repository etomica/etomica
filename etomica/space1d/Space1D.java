package etomica.space1d;

import java.io.ObjectStreamException;

import etomica.EtomicaInfo;
import etomica.NearestImageTransformer;
import etomica.Space;

public final class Space1D extends Space {
    
    /**
     * Private constructor for singleton.
     */
    private Space1D() {
        super();
    }
    
    /**
     * @return Returns the instance.
     */
    public static Space1D getInstance() {
        return INSTANCE;
    }

    public final int D() {
        return 1;
    }

    public final int powerD(int n) {
        return n;
    }

    public final double powerD(double a) {
        return a;
    }
    
    /**
     * Returns the Dth root of the given value, a^(1/D), where D is the dimension of the space.
     */
    public double rootD(double a) {return a;}

    public double sphereVolume(double r) {
        return 2.0 * r;
    } //volume of a sphere of radius r

    public double sphereArea(double r) {
        return 2.0;
    } //surface area of sphere of radius r (used for differential shell volume)

    public etomica.space.Vector makeVector() {
        return new Vector1D();
    }

    public etomica.space.Orientation makeOrientation() {
        System.out.println("Orientation class not implemented in 1D");
        return null;
    }

    public etomica.space.Tensor makeTensor() {
        return new Tensor1D();
    }

    public etomica.space.Tensor makeRotationTensor() {
        return new RotationTensor();
    }

    public int[] makeArrayD(int i) {
        return new int[] { i };
    }

    public double[] makeArrayD(double d) {
        return new double[] { d };
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("One-dimensional space");
        return info;
    }

    /**
     * Computes the square of the magnitude of the difference of two vectors, subject
     * to a nearest image transformation.  This method constructs a new vector that
     * is used as the work-vector input to the other r2 method.
     */
    public static final double r2(Vector1D u1, Vector1D u2, NearestImageTransformer b) {
        return r2(u1, u2, b, new Vector1D());
    }

    /**
     * Computes the square of the magnitude of the difference of two vectors, subject
     * to a nearest image transformation.  
     * @param u1 one of the vectors and is unchanged by the calculation
     * @param u2 the other vector and is unchanged by the calculation
     * @param b a nearest image transformation
     * @param work a work vector used for the calculation.
     */
    public static final double r2(Vector1D u1, Vector1D u2, NearestImageTransformer b,
            Vector1D work) {
        work.Ev1Mv2(u1, u2);
        b.nearestImage(work);
        return work.squared();
    }
    
    /**
     * Required to guarantee singleton when deserializing.
     * @return the singleton INSTANCE
     */
    private Object readResolve() throws ObjectStreamException {
        return INSTANCE;
    }
    
    private static final Space1D INSTANCE = new Space1D();

}
