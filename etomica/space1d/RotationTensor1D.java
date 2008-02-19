package etomica.space1d;

public class RotationTensor1D extends Tensor1D implements etomica.space.RotationTensor {

    private static final long serialVersionUID = 1L;
    public RotationTensor1D() {super(); reset();}
    public void reset() {
        xx = 1.0;
    }
    public void setAxial(int i, double theta) {
    }
    public void setAngles(double[] angles) {}
}
