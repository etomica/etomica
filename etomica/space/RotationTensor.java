package etomica.space;

/**
 * A tensor intended for use to define rotations in space.  Includes methods
 * to set its components in terms of rotation angles.
 */
public interface RotationTensor extends Tensor {
    
    /**
     * Sets tensor to a condition of no rotation.
     */
    public void reset();

    public void invert();

    public void setAxial(int i, double theta);

    public void setAngles(double[] angles);
}