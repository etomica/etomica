/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.math.function.IFunction;

/**
 * Interface for a second-order tensor. Implementations differ depending on the dimension of
 * the space where the tensor is used.
 */
public interface Tensor extends Cloneable {

    /**
     * @return a new tensor that is a copy of this one
     */
    Object clone();

    /**
     * @return the dimension of the space in which this tensor is defined
     */
    int D();

    /**
     * @param i "row" index of the requested component
     * @param j "column" index of the requested component
     * @return the component a<sub>i,j</sub> of this tensor
     */
    double component(int i, int j);

    /**
     * @param i "row" index of the component being set
     * @param j "column" index of the component being set
     * @param d value to which the component a<sub>i,j</sub> of this tensor is set
     */
    void setComponent(int i, int j, double d);

    /**
     * "Equals" operation. Sets all components of this tensor equal to those of the given tensor.
     * @param t the given tensor used to set this tensor's values
     */
    void E(Tensor t);
    
    /**
     * Fills the tensor column-wise with the given vectors.  The number
     * of vectors must equal the dimension of the tensor, and the vector
     * dimensions must equal the tensor dimension. Each vector forms a column of the tensor, so
     * for D = 3, xx = v[0].x, xy = v[1].x, xz = v[2].x,...,yx = v[0].y, etc.
     * @param v the given vectors used to set this tensor's values
     */
    void E(Vector[] v);

    /**
     * Sets the tensor elements using the elements of the array.
     * Tensor values are filled row-wise, so array values [0][0], [0][1], [0][2],...[1][0]... are
     * assigned to xx, xy, xz, yx, etc. respectively.
     */
    void E(double[][] d);

    /**
     * "Equals" operation, setting all elements of this tensor to the given value
     * @param a the given value used to set this tensor's elements
     */
    void E(double a);

    /**
     * Makes this a diagonal tensor with elements set by the given vector. All
     * other elements of this tensor are set to zero.
     * @param v the given vector
     */
    void diagE(Vector v);

    /**
     * Makes this a diagonal tensor with elements set all equal to the given value. All
     * other elements of this tensor are set to zero.
     * @param a the given value
     */
    void diagE(double a);

    /**
     * Assigns the tensor elements column-wise to the given vectors. The number
     * of vectors must equal the dimension of the tensor, and the vector
     * dimensions must equal the tensor dimension.
     */
    void assignTo(Vector[] v);
    
    /**
     * Sets this tensor equal to the dyadic or outer product of the given vectors.
     * Element ab is given by v1.a * v2.b
     * @param v1 first vector of the dyad
     * @param v2 second vector of the dyad
     */
    void Ev1v2(Vector v1, Vector v2);

    /**
     * "Plus equals" (+=) operation, adding the given value to each element of this tensor
     * @param a the value added to each element
     */
    void PE(double a);

    /**
     * "Plus equals" (+=) operation, adding the elements of the given tensor to the
     * respective elements of this tensor
     * @param t the given tensor added to this tensor
     */
    void PE(Tensor t);

    /**
     * "Plus equals" (+=) operation, performed on a single element of this tensor.
     *
     * @param i "row" index of the modified element
     * @param j "column" index of the modified element
     * @param a value to be added to the element a<sub>i,j</sub>
     */
    void PE(int i, int j, double a);
    
    /**
     * "Plus equals" (+=) operation with a dyad argument. Increments this tensor by the dyadic or outer product of the given vectors.
     * Element ab of this tensor is replaced by (ab + v1.a * v2.b)
     */
    void PEv1v2(Vector v1, Vector v2);


    /**
     * "Minus equals" (-=) operation with a dyad argument. Decrements this tensor by the dyadic or outer product of the given vectors.
     * Element ab of this tensor is replaced by (ab - v1.a * v2.b)
     */
    void MEv1v2(Vector v1, Vector v2);


    /**
     * "Minus equals" (-=) operation, subtracting the elements of the given tensor from the
     * respective elements of this tensor. This tensor's values are replaced by the so-computed difference.
     * @param t the given tensor subtracted from this tensor
     */
    void ME(Tensor t);

    /**
     * "Plus equals" (+=) operation based on a scalar and tensor. Each element ab of this tensor is
     * replaced by a1 * t1.ab.
     * @param a1 the scalar multiplier
     * @param t1 the given tensor
     */
    void PEa1Tt1(double a1, Tensor t1);

    /**
     *
     * @return the sum of the diagonal elements of this tensor
     */
    double trace();

    /**
     * Performs a transpose operation on this tensor, replacing each element ab by its counterpart element ba
     */
    void transpose();

    /**
     * Computes the determinant of the matrix represented by this tensor
     * @return the determinant of this tensor
     */
    double determinant();

    /**
     * Replaces this tensor with its inverse. The inverse is defined such that the tensor multiplied by its inverse
     * gives the identity tensor.
     */
    void invert();


    /**
     * "Times equals" (*=) operation, replaces each element of this tensor by its current value times the given scalar
     * @param a the given scalar multiplier
     */
    void TE(double a);

    /**
     * "Times equals" (*=) operation, as matrix multiplication. Replaces each element of this tensor by the
     * value given by result of (this)*t, where multiplication is treated as if multiplying matrices.
     * @param t the tensor this is multiplying
     */
    void TE(Tensor t);

    /**
     * "Divide equals" (/=) operation, element-by-element. Replaces each element of this tensor by its current value divided by the
     * corresponding element in the given tensor. Element ab is replated by ab / t.ab.
     * @param t the given tensor
     */
    void DE(Tensor t);

    /**
     * Returns tensor in a new array, writing elements row-wise. For example, in 2D: xx, xy, yx, yy.
     * @return a new array containing the tensor elements
     */
    double[] toArray();

    /**
     * @return true if any element of this tensor is NaN (not a number)
     */
    boolean isNaN();

    /**
     * Applies given function to each element of this tensor. Replaces element with the resulting value.
     * @param f the function applied to each element
     */
    void map(IFunction f);

    /**
     * Applies the given tensor transformation to this vector, replacing its
     * elements with the transformed values. This operation is equivalent to a matrix
     * multiplication.
     */
    void transform(Vector x);

    /**
     * Fills the given array with the elements of this tensor.
     * Array values 0, 1, 2,... are assigned from xx, xy, xz, yx, etc. respectively,
     * filling in a row-wise fashion.
     * @param d an array of length D<sup>2</sup>, containing the tensor elements upon return
     */
    void assignTo(double[] d);

    /**
     * Fills the given 2-dimensional array with the elements of this tensor.
     * Array values [0][0], [0][1] are assigned from xx, xy, xz, yx, etc. respectively,
     * such that d[i] contains the i<sup>th</sup> row of the tensor.
     * @param d a D-by-D array, containing the tensor elements upon return
     */
    void assignTo(double[][] d);

    /**
     * Element-by-element comparison of tensors
     * @return true if every element of given tensor equals the corresponding element of this tensor
     * @param t the given tensor
     */
    boolean equals(Tensor t);
}
