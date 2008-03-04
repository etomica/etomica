package etomica.space;

import etomica.api.IFunction;
import etomica.api.IVector;


public interface Tensor extends Cloneable {
    public Object clone();
    
    public abstract int D();
    public abstract double component(int i, int j);
    public abstract void setComponent(int i, int j, double d);
    public abstract void E(Tensor t);
    
    /**
     * Fills the tensor column-wise with the given vectors.  The number
     * of vectors must equal the dimension of the tensor, and the vector
     * dimensions must equal the tensor dimension.
     * @param v
     */
    public abstract void E(IVector[] v);
    
    /**
     * Assigns the tensor elements column-wise to the given vectors. The number
     * of vectors must equal the dimension of the tensor, and the vector
     * dimensions must equal the tensor dimension.
     */
    public abstract void assignTo(IVector[] v);
    
    /**
     * Sets this equal to the dyadic or outer product of the given vectors.
     * Element ab is given by v1.a * v2.b 
     */
    public abstract void Ev1v2(IVector v1, IVector v2);
    
    public abstract void E(double a);
    public abstract void PE(double a);
    public abstract void PE(Tensor t);
    public abstract void PE(int i, int j, double a);
    
    /**
     * Increments this by the dyadic or outer product of the given vectors.
     */
    public abstract void PEv1v2(IVector v1, IVector v2);
    public abstract void MEv1v2(IVector v1, IVector v2);
    public abstract void ME(Tensor t);
    public abstract void PEa1Tt1(double a1, Tensor t1);
    public abstract double trace();
    public abstract void transpose();
    public abstract void invert();
    public abstract void TE(double a);
    public abstract void TE(Tensor t);
    public abstract void DE(Tensor t);
    public double[] toArray();
    public boolean isNaN();
    public void map(IFunction f);
    /**
     * Sets the tensor elements using the elements of the array.
     * Tensor values are filled row-wise, so array values 0, 1, 2,... are 
     * assigned to xx, xy, xz, yx, etc. respectively.
     */
    public abstract void E(double[] d);
    
    /**
     * Applies the given tensor transformation to this vector, replaced its
     * elements with the transformed values.
     */
    public void transform(IVector A);

    /**
     * Fills the given array with the elements of this tensor.
     * Array values 0, 1, 2,... are assigned from xx, xy, xz, yx, etc. respectively,
     * filling in a row-wise fashion.
     * @param d
     */
    public abstract void assignTo(double[] d);
    
    public abstract void assignTo(double[][] d);
}