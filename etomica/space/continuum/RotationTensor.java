package etomica.space.continuum;
import etomica.Space;

public class RotationTensor extends Tensor implements Space.RotationTensor {
    public RotationTensor() {super(); reset();}
    public void reset() {
        xx = 1.0; xy = 0.0; xz = 0.0;
        yx = 0.0; yy = 1.0; yz = 0.0;
        zx = 0.0; zy = 0.0; zz = 1.0;
    }
    public void setAxial(int i, double theta) {
        double st = Math.sin(theta);
        double ct = Math.cos(theta);
        switch(i) {
            case 0: xx = 1.; xy = 0.; xz = 0.;
                    yx = 0.; yy = ct; yz = -st;
                    zx = 0.; zy = st; zz = ct;
                    return;
            case 1: xx = ct; xy = 0.; xz = -st;
                    yx = 0.; yy = 1.; yz = 0.;
                    zx = st; zy = 0.; zz = ct;
                    return;
            case 2: xx = ct; xy = -st; xz = 0.;
                    yx = st; yy = ct;  yz = 0.;
                    zx = 0.; zy = 0.;  zz = 1.;
                    return;
            default: throw new IllegalArgumentException();
        }
    }
    public void setAngles(double[] angles) {}
    public void invert() {}
}

