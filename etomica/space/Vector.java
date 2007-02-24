package etomica.space;

import etomica.util.Function;

/**
 * Class defining methods for a vector in a D-dimensional space.  Arithmetic vector operations
 * are named with capital letters that are coded to the operations performed:
 * <ul>
 * <li>E for "equals" =
 * <li>P for "plus", +
 * <li>M for "minus", -
 * <li>T for "times", *
 * <li>D for "divided", /
 * </ul>
 * The operation may be interpreted by writing the corresponding symbolic expression; for
 * example, PE is "+=" and will cause each element of this vector to be incremented by
 * some amount (specific to the method).  Most arithemetic operations involving two
 * or more vectors simply apply element-by-element.  Additionally, standard vector
 * operations are defined, such as the dot product and cross product (indicatate with "X" in the name).
 * <p>
 * Vector names also may include a1, a2, v1, v2, etc.  An "a" indicates a scalar, while
 * a "v" indicates another vector.  Thus the method Ea1Tv1 (read "= a1*v1) will replace
 * all elements of this vector with the elements of vector v1, each multiplied by the
 * scalar a1.
 * <p>
 * Vector components may be accessed and changed via the methods "x" and "setX", respectively. 
 * 
 */
public abstract class Vector implements java.io.Serializable {

    /**
     * Dimension of the space occupied by the vector. Number of elements in the
     * vector.
     */
    public abstract int D();

    /**
     * Assigns the components of this vector to the elements of the given array.
     * Does not check that array length is consistent with vector dimension.
     * Inverse of the E method.
     */
    public abstract void assignTo(double[] array);

    /**
     * Returns a new array with elements equal to the components of this vector.
     */
    public abstract double[] toArray();

    /**
     * Returns true if all corresponding elements of this and the given vector
     * are equal; returns false otherwise.
     */
    public abstract boolean equals(Vector v);

    /**
     * Computes the spherical-coordinate representation of this vector, and
     * returns them in the given array.
     */
    public abstract void sphericalCoordinates(double[] result);

    /**
     * Vector components corresponding to the given index. For example,
     * x-component is given for i = 0. If index exceeds dimension of vector, no
     * error is given and the last element of the vector is returned.
     */
    public abstract double x(int i);

    /**
     * Sets the i-th component of the vector to the given value d. If index
     * exceeds dimension of vector, no error is given and the last element of
     * the vector is assigned.
     */
    public abstract void setX(int i, double d); 
    
    /**
     * Sets the components of this vector equal to those of the given vector.
     */
    public abstract void E(Vector u);

    /**
     * Sets all components of this vector equal to the given value.
     */
    public abstract void E(double a);

    /**
     * Sets the components of this vector equal to the corresponding elements of
     * the given array. Inverse of the assignTo method.
     */
    public abstract void E(double[] a);

    /**
     * Sets the components of this vector equal to the corresponding elements of
     * the given integer array.
     */
    public abstract void E(int[] a);

    /**
     * Plus-equals (+=) operation. Increments each component of this vector by
     * the corresponding value in the given vector.
     */
    public abstract void PE(Vector u);

    /**
     * Plus-equals (+=) operation applied to one component.
     * 
     * @param i
     *            index of the component being incremented
     * @param a
     *            increment value
     */
    public abstract void PE(int i, double a);

    /**
     * Plus-equals (+=) operation, causing a constant value to be added to all
     * components of this vector.
     */
    public abstract void PE(double a); //adds a constant value to all elements

    /**
     * Minus-equals (-=) operation. Decrements each component of this vector by
     * the corresponding value in the given vector.
     */
    public abstract void ME(Vector u);

    /**
     * Times-equals (*=) operation. Multiplies each component of this vector by
     * the corresponding value in the given vector.
     */
    public abstract void TE(Vector u);

    /**
     * Divide-equals (/=) operation. Divides each component of this vector by
     * the corresponding value in the given vector.
     */
    public abstract void DE(Vector u);

    /**
     * Multiplies all components by a constant.
     */
    public abstract void TE(double a); //multipies all components by the

    // constant a

    /**
     * Multiplies one component, of index i, by a constant.
     */
    public abstract void TE(int i, double a);

    /**
     * Operation (= a1 * v1); sets the components of this to those of the given
     * vector multiplied by the given constant.
     */
    public abstract void Ea1Tv1(double a, Vector v1); //sets this vector to a*u

    /**
     * Sets this vector equal to (v1 + a1 * v2).
     */
    public abstract void Ev1Pa1Tv2(Vector v1, double a1, Vector v2);

    /**
     * Increments (+=) components by (a1 * v1)
     */
    public abstract void PEa1Tv1(double a, Vector v1);

    /**
     * Sets the components of this vector equal to (v1 + v2)
     */
    public abstract void Ev1Pv2(Vector v1, Vector v2);

    /**
     * Adds/subtracts a scalar to each element, with the choice of 
     * addition or subtraction for each element determined by the corresponding
     * element of the vector v1.  In symbols, this is += a1*sgn(v1).  Zero
     * is added if the corresponding element of v1 is zero.
     */
    public abstract void PEa1SGNv1(double a1, Vector v1);
    
    // vectors
    /**
     * Sets the components of this vector equal to (v1 - v2)
     */
    public abstract void Ev1Mv2(Vector v1, Vector v2);

    /**
     * Adds (Plus, +) this vector to the given one and returns sum in a new vector.
     */
    public abstract Vector P(Vector u);

    /**
     * Subtracts (Minus, -) given vector from this returns result in a new vector.
     */
    public abstract Vector M(Vector u);

    /**
     * Multiplies (Times, *) given vector times this and returns result 
     * in a new vector.
     */
    public abstract Vector T(Vector u); 

    /**
     * Divides (/) this vector by the given one and returns result in a new vector.
     */
    public abstract Vector D(Vector u);

    /**
     * Returns square of vector resulting from subtracting given vector
     * from this one.  Neither vector is changed by this operation.
     */
    public abstract double Mv1Squared(Vector v1);

    /**
     * Sets to zero all elements having magnitude (absolute value) less than the given value.
     */
    public abstract void truncate(double eps);
    
    /**
     * Replaces each component of this vector with its absolute value.
     */
    public abstract void abs();

    /**
     * Replaces each component of this vector with its value modulo a
     */
    public abstract void mod(double a); 

    /**
     * Replaces each component of this vector with its value modulo
     * the corresponding component in u
     */
    public abstract void mod(Vector u); //each component replaced with itself

    /**
     * Sets this equal to (r mod u) - r.  Fails if this == r.
     */
    public abstract void EModShift(Vector r, Vector u);

    /**
     * Hard to explain.
     */
    public abstract void EMod2Shift(Vector r, Vector u);

    /**
     * Returns the minimum of all components of this vector.
     */
    public abstract double min();

    /**
     * Returns the maximum of all components of this vector.
     */
    public abstract double max();

    /**
     * Replaces each element of this vector with the result of taking the minimum of it and the
     * corresponding element of the given vector.  If the corresponding element is less
     * (closer to negative infinity), it replace the element of this vector.  Operation
     * may change all, some, or none of this vector's elements.  
     */
    public abstract void minE(Vector v);
    
    /**
     * Replaces each element of this vector with the result of taking the maximum of it and the
     * corresponding element of the given vector.  If the corresponding element is greater
     * (closer to positive infinity), it replace the element of this vector.  Operation 
     * may change all, some, or none of this vector's elements.  
     */
    public abstract void maxE(Vector v);
    
    /**
     * Returns the square magnitude of this vector, e.g., x^2 + y^2 for D = 2.
     */
    public abstract double squared();

    /**
     * Normalizes this vector, so this.squared() == 1.  Divides all 
     * components by Math.sqrt(this.squared()).
     */
    public abstract void normalize();

    /**
     * Returns true if all components of this vector are zero; false otherwise.
     */
    public abstract boolean isZero();

    /**
     * Returns the dot product of this vector with the given one.
     */
    public abstract double dot(Vector u); 

    /**
     * Applies the given tensor transformation to this vector, replaced its
     * elements with the transformed values.
     */
    public abstract void transform(Tensor A);

    /**
     * Assigns each components to (its own) random value between 0 and d.
     */
    public abstract void setRandom(double d);

    /**
     * Assigns each component to (its own) random value between 0 and v(i)
     */
    public abstract void setRandom(Vector v);
    
    /**
     * Assigns this vector to equal a point chosen randomly on the 
     * surface of a unit sphere.
     */
    public abstract void setRandomSphere();

    /**
     * Assigns each component to (its own) random value between -0.5 and + 0.5.
     */
    public abstract void setRandomCube();

    /**
     * Assigns this vector to equal a point chosen randomly in the volume
     * of a unit spheres.
     */
    public abstract void setRandomInSphere();

    /**
     * Assigns this vector to a new one chosen randomly in the cube of
     * size 2*d (selecting +/- d) centered on its present position.
     */
    public abstract void randomStep(double d); 
    
    /**
     * Performs a random rotation of this vector within a cone of 
     * angle theta on its current position.
     */
    public abstract void randomRotate(double thetaStep);
    
    /**
     * Returns true if any element of the vector is not-a-number.
     */
    public abstract boolean isNaN();
    
    /**
     * Applies the given function to each element of the vector.
     */
    public abstract void map(Function f);
}