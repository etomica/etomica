package etomica.space1d;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class RotationTensor1D extends Tensor1D implements etomica.space.RotationTensor {
    public RotationTensor1D() {super(); reset();}
    public void reset() {
        xx = 1.0;
    }
    public void setAxial(int i, double theta) {
    }
    public void setAngles(double[] angles) {}
    public void invert() {}
}
