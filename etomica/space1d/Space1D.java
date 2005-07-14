package etomica.space1d;

import java.io.ObjectStreamException;

import etomica.EtomicaInfo;
import etomica.Space;
import etomica.space.Boundary;

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

    public static final double r2(Vector1D u1, Vector1D u2, Boundary b) {
        Vector1D.WORK.x = u1.x - u2.x;
        b.nearestImage(Vector1D.WORK);
        return Vector1D.WORK.x * Vector1D.WORK.x;
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
