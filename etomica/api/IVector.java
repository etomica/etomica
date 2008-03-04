package etomica.api;


/**
 * Basic vector interface containing methods needed for basic vector operations
 * and other commonly used methods.
 */
public interface IVector {

    /**
     * Dimension of the space occupied by the vector. Number of elements in the
     * vector.
     */
    public int getD();

    /**
     * Assigns the components of this vector to the elements of the given array.
     * Does not check that array length is consistent with vector dimension.
     * Inverse of the E method.
     */
    public void assignTo(double[] array);

    /**
     * Returns true if all corresponding elements of this and the given vector
     * are equal; returns false otherwise.
     */
    public boolean equals(IVector v);

    /**
     * Vector components corresponding to the given index. For example,
     * x-component is given for i = 0. If index exceeds dimension of vector, no
     * error is given and the last element of the vector is returned.
     */
    public double x(int i);

    /**
     * Sets the i-th component of the vector to the given value d. If index
     * exceeds dimension of vector, no error is given and the last element of
     * the vector is assigned.
     */
    public void setX(int i, double d);

    /**
     * Sets the components of this vector equal to those of the given vector.
     */
    public void E(IVector u);

    /**
     * Sets all components of this vector equal to the given value.
     */
    public void E(double a);

    /**
     * Sets the components of this vector equal to the corresponding elements of
     * the given array. Inverse of the assignTo method.
     */
    public void E(double[] a);

    /**
     * Plus-equals (+=) operation. Increments each component of this vector by
     * the corresponding value in the given vector.
     */
    public void PE(IVector u);

    /**
     * Plus-equals (+=) operation, causing a constant value to be added to all
     * components of this vector.
     */
    public void PE(double a); //adds a constant value to all elements

    /**
     * Minus-equals (-=) operation. Decrements each component of this vector by
     * the corresponding value in the given vector.
     */
    public void ME(IVector u);

    /**
     * Times-equals (*=) operation. Multiplies each component of this vector by
     * the corresponding value in the given vector.
     */
    public void TE(IVector u);

    /**
     * Divide-equals (/=) operation. Divides each component of this vector by
     * the corresponding value in the given vector.
     */
    public void DE(IVector u);

    /**
     * Multiplies all components by a constant.
     */
    public void TE(double a);

    /**
     * Operation (= a1 * v1); sets the components of this to those of the given
     * vector multiplied by the given constant.
     */
    public void Ea1Tv1(double a, IVector v1);

    /**
     * Increments (+=) components by (a1 * v1)
     */
    public void PEa1Tv1(double a, IVector v1);

    /**
     * Sets the components of this vector equal to (v1 + v2)
     */
    public void Ev1Pv2(IVector v1, IVector v2);

    /**
     * Sets the components of this vector equal to (v1 - v2)
     */
    public void Ev1Mv2(IVector v1, IVector v2);

    /**
     * Replaces each component of this vector with its value modulo
     * the corresponding component in u
     */
    public void mod(IVector u);

    /**
     * Returns the square magnitude of this vector, e.g., x^2 + y^2 for D = 2.
     */
    public double squared();

    /**
     * Returns true if all components of this vector are zero; false otherwise.
     */
    public boolean isZero();

    /**
     * Returns the dot product of this vector with the given one.
     */
    public double dot(IVector u);

    /**
     * Returns true if any element of the vector is not-a-number.
     */
    public boolean isNaN();

    /**
     * Returns square of vector resulting from subtracting given vector
     * from this one.  Neither vector is changed by this operation.
     */
    public double Mv1Squared(IVector v1);

    /**
     * Normalizes this vector, so this.squared() == 1.  Divides all 
     * components by Math.sqrt(this.squared()).
     */
    public void normalize();

    /**
     * Applies the given function to each element of the vector.
     */
    public void map(IFunction f);

}