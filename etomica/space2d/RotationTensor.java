package etomica.space2d;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class RotationTensor extends Tensor2D implements etomica.space.RotationTensor {
    public RotationTensor() {super(); reset();}
    public void reset() {
        xx = 1.0; xy = 0.0;
        yx = 0.0; yy = 1.0;
    }
    public void setAxial(int i, double theta) {
        double st = Math.sin(theta);
        double ct = Math.cos(theta);
        switch(i) {
            case 2: xx = ct; xy=-st;
                    yx = st; yy=ct;
                    return;
            default: throw new IllegalArgumentException();
        }
    }
    public void setAngles(double[] angles) {}
    public void invert() {xy *= -1; yx *= -1;}
}
