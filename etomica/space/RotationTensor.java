package etomica.space;

/**
 * A tensor intended for use to define rotations in space.  Includes methods
 * to sets its components in terms of rotation angles.
 */
public interface RotationTensor extends Tensor {
    
    public void reset();

    public void invert();

    public void setAxial(int i, double theta);

    public void setAngles(double[] angles);
}