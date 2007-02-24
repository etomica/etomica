package etomica.space;

import etomica.util.Function;

public interface IVector {

    /**
     * Dimension of the space occupied by the vector. Number of elements in the
     * vector.
     */
    public int D();

    /**
     * Assigns the components of this vector to the elements of the given array.
     * Does not check that array length is consistent with vector dimension.
     * Inverse of the E method.
     */
    public void assignTo(double[] array);

    /**
     * Returns a new array with elements equal to the components of this vector.
     */
    public double[] toArray();

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
     * Sets the components of this vector equal to the corresponding elements of
     * the given integer array.
     */
    public void E(int[] a);

    /**
     * Plus-equals (+=) operation. Increments each component of this vector by
     * the corresponding value in the given vector.
     */
    public void PE(IVector u);

    /**
     * Plus-equals (+=) operation applied to one component.
     * 
     * @param i
     *            index of the component being incremented
     * @param a
     *            increment value
     */
    public void PE(int i, double a);

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
    public void TE(double a); //multipies all components by the

    /**
     * Multiplies one component, of index i, by a constant.
     */
    public void TE(int i, double a);

    /**
     * Operation (= a1 * v1); sets the components of this to those of the given
     * vector multiplied by the given constant.
     */
    public void Ea1Tv1(double a, IVector v1); //sets this vector to a*u

    /**
     * Sets this vector equal to (v1 + a1 * v2).
     */
    public void Ev1Pa1Tv2(IVector v1, double a1, IVector v2);

    /**
     * Increments (+=) components by (a1 * v1)
     */
    public void PEa1Tv1(double a, IVector v1);

    /**
     * Sets the components of this vector equal to (v1 + v2)
     */
    public void Ev1Pv2(IVector v1, IVector v2);

    /**
     * Adds/subtracts a scalar to each element, with the choice of 
     * addition or subtraction for each element determined by the corresponding
     * element of the vector v1.  In symbols, this is += a1*sgn(v1).  Zero
     * is added if the corresponding element of v1 is zero.
     */
    public void PEa1SGNv1(double a1, IVector v1);

    // vectors
    /**
     * Sets the components of this vector equal to (v1 - v2)
     */
    public void Ev1Mv2(IVector v1, IVector v2);

    /**
     * Adds (Plus, +) this vector to the given one and returns sum in a new vector.
     */
    public IVector P(IVector u);

    /**
     * Subtracts (Minus, -) given vector from this returns result in a new vector.
     */
    public IVector M(IVector u);

    /**
     * Multiplies (Times, *) given vector times this and returns result 
     * in a new vector.
     */
    public IVector T(IVector u);

    /**
     * Divides (/) this vector by the given one and returns result in a new vector.
     */
    public IVector D(IVector u);

    /**
     * Returns square of vector resulting from subtracting given vector
     * from this one.  Neither vector is changed by this operation.
     */
    public double Mv1Squared(IVector v1);

    /**
     * Sets to zero all elements having magnitude (absolute value) less than the given value.
     */
    public void truncate(double eps);

    /**
     * Replaces each component of this vector with its absolute value.
     */
    public void abs();

    /**
     * Replaces each component of this vector with its value modulo a
     */
    public void mod(double a);

    /**
     * Replaces each component of this vector with its value modulo
     * the corresponding component in u
     */
    public void mod(IVector u); //each component replaced with itself

    /**
     * Sets this equal to (r mod u) - r.  Fails if this == r.
     */
    public void EModShift(IVector r, IVector u);

    /**
     * Hard to explain.
     */
    public void EMod2Shift(IVector r, IVector u);

    /**
     * Returns the minimum of all components of this vector.
     */
    public double min();

    /**
     * Returns the maximum of all components of this vector.
     */
    public double max();

    /**
     * Replaces each element of this vector with the result of taking the minimum of it and the
     * corresponding element of the given vector.  If the corresponding element is less
     * (closer to negative infinity), it replace the element of this vector.  Operation
     * may change all, some, or none of this vector's elements.  
     */
    public void minE(IVector v);

    /**
     * Replaces each element of this vector with the result of taking the maximum of it and the
     * corresponding element of the given vector.  If the corresponding element is greater
     * (closer to positive infinity), it replace the element of this vector.  Operation 
     * may change all, some, or none of this vector's elements.  
     */
    public void maxE(IVector v);

    /**
     * Returns the square magnitude of this vector, e.g., x^2 + y^2 for D = 2.
     */
    public double squared();

    /**
     * Normalizes this vector, so this.squared() == 1.  Divides all 
     * components by Math.sqrt(this.squared()).
     */
    public void normalize();

    /**
     * Returns true if all components of this vector are zero; false otherwise.
     */
    public boolean isZero();

    /**
     * Returns the dot product of this vector with the given one.
     */
    public double dot(IVector u);

    /**
     * Applies the given tensor transformation to this vector, replaced its
     * elements with the transformed values.
     */
    public void transform(Tensor A);

    /**
     * Computes the spherical-coordinate representation of this vector, and
     * returns them in the given array.
     */
    public void sphericalCoordinates(double[] result);

    /**
     * Assigns each components to (its own) random value between 0 and d.
     */
    public void setRandom(double d);

    /**
     * Assigns each component to (its own) random value between 0 and v(i)
     */
    public void setRandom(IVector v);

    /**
     * Assigns this vector to equal a point chosen randomly on the 
     * surface of a unit sphere.
     */
    public void setRandomSphere();

    /**
     * Assigns each component to (its own) random value between -0.5 and + 0.5.
     */
    public void setRandomCube();

    /**
     * Assigns this vector to equal a point chosen randomly in the volume
     * of a unit spheres.
     */
    public void setRandomInSphere();

    /**
     * Assigns this vector to a new one chosen randomly in the cube of
     * size 2*d (selecting +/- d) centered on its present position.
     */
    public void randomStep(double d);

    /**
     * Performs a random rotation of this vector within a cone of 
     * angle theta on its current position.
     */
    public void randomRotate(double thetaStep);

    /**
     * Returns true if any element of the vector is not-a-number.
     */
    public boolean isNaN();

    /**
     * Applies the given function to each element of the vector.
     */
    public void map(Function f);

}