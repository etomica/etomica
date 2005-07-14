package etomica.space3d;

import etomica.EtomicaInfo;
import etomica.Space;
import etomica.space.Boundary;

/**
 * 
 * Factory and methods appropriate to a 3-dimensional space.
 * 
 * @author David Kofke
 *
 */
public class Space3D extends Space {

    public final int D() {
        return D;
    }

    public final int powerD(int n) {
        return n * n * n;
    }

    public final double powerD(double a) {
        return a * a * a;
    }

    public int[] makeArrayD(int i) {
        return new int[] { i, i, i };
    }

    public double[] makeArrayD(double d) {
        return new double[] { d, d, d };
    }

    public static final Space3D INSTANCE = new Space3D();

    public Space3D() {
        super(3);
    }

    public double sphereVolume(double r) {
        return (Math.PI * 4.0 * r * r * r / 3.0);
    }

    public double sphereArea(double r) {
        return (Math.PI * 4 * r * r);
    }

    public etomica.space.Vector makeVector() {
        return new Vector3D();
    }

    public etomica.space.Orientation makeOrientation() {
        return new Orientation();
    }

    public etomica.space.Tensor makeTensor() {
        return new Tensor3D();
    }

    public etomica.space.Tensor makeRotationTensor() {
        return new RotationTensor();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Three-dimensional space");
        return info;
    }

    public static final double r2(Vector3D u1, Vector3D u2, Boundary b) {
        return r2(u1, u2, b, new Vector3D());
    }

    public static final double r2(Vector3D u1, Vector3D u2, Boundary b,
            Vector3D work) {
        work.Ev1Mv2(u1, u2);
        b.nearestImage(work);
        return work.squared();
    }

}//end of Space3D
