package etomica.space;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public interface RotationTensor extends Tensor {
    public void reset();
    public void invert();
    public void setAxial(int i, double theta);
    public void setAngles(double[] angles);
}