package etomica.space;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public interface Tensor {
    public abstract int length();
    public abstract double component(int i, int j);
    public abstract void setComponent(int i, int j, double d);
    public abstract void E(Tensor t);
    public abstract void E(Vector u1, Vector u2);
    public abstract void E(double a);
    public abstract void PE(Tensor t);
    public abstract void PE(int i, int j, double a);
    public abstract void PE(Vector u1, Vector u2);
    public abstract double trace();
    public abstract void TE(double a);
    /**
     * Sets the tensor elements using the elements of the array.
     * Array values 0, 1, 2,... are assigned to xx, xy, xz, yx, etc. respectively.
     */
    public abstract void E(double[] d);
    
    /**
     * Fills the given array with the elements of this tensor.
     * Array values 0, 1, 2,... are assigned from xx, xy, xz, yx, etc. respectively.
     * @param d
     */
    public abstract void assignTo(double[] d);
}