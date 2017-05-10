package etomica.space1d;

import etomica.space.RotationTensor;

/**
 * Rotation tensor in 1D either flips orientation or keeps it unchanged.
 * @author David Kofke
 */
public class RotationTensor1D extends Tensor1D implements RotationTensor {

    /**
     * Default is identity tensor that performs no rotation (flip).
     */
    public RotationTensor1D() {
        super(1.0);
    }

    public void reset() {
        xx = 1.0;
    }

    /**
     * @param i ignored
     * @param theta sets tensor to do flip if nonzero; no flip if theta is zero
     */
    public void setAxial(int i, double theta) {
        xx = (theta == 0.0) ? +1 : -1;
    }
}
