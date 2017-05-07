/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.api.IFunction;
import etomica.api.IVector;
import etomica.api.IVectorMutable;


public interface Tensor extends Cloneable {
    Object clone();
    
    int D();
    double component(int i, int j);
    void setComponent(int i, int j, double d);
    void E(Tensor t);
    
    /**
     * Fills the tensor column-wise with the given vectors.  The number
     * of vectors must equal the dimension of the tensor, and the vector
     * dimensions must equal the tensor dimension.
     * @param v
     */
    void E(IVector[] v);
    
    /**
     * Assigns the tensor elements column-wise to the given vectors. The number
     * of vectors must equal the dimension of the tensor, and the vector
     * dimensions must equal the tensor dimension.
     */
    void assignTo(IVectorMutable[] v);
    
    /**
     * Sets this equal to the dyadic or outer product of the given vectors.
     * Element ab is given by v1.a * v2.b 
     */
    void Ev1v2(IVector v1, IVector v2);
    
    void E(double a);
    void PE(double a);
    void PE(Tensor t);
    void PE(int i, int j, double a);
    
    /**
     * Increments this by the dyadic or outer product of the given vectors.
     */
    void PEv1v2(IVector v1, IVector v2);
    void MEv1v2(IVector v1, IVector v2);
    void ME(Tensor t);
    void PEa1Tt1(double a1, Tensor t1);
    double trace();
    void transpose();
    double determinant();
    void invert();
    void TE(double a);
    void TE(Tensor t);
    void DE(Tensor t);
    double[] toArray();
    boolean isNaN();
    void map(IFunction f);
    /**
     * Sets the tensor elements using the elements of the array.
     * Tensor values are filled row-wise, so array values 0, 1, 2,... are 
     * assigned to xx, xy, xz, yx, etc. respectively.
     */
    void E(double[][] d);

    /**
     * Applies the given tensor transformation to this vector, replaced its
     * elements with the transformed values.
     */
    void transform(IVectorMutable A);

    /**
     * Fills the given array with the elements of this tensor.
     * Array values 0, 1, 2,... are assigned from xx, xy, xz, yx, etc. respectively,
     * filling in a row-wise fashion.
     * @param d
     */
    void assignTo(double[] d);
    
    void assignTo(double[][] d);
}