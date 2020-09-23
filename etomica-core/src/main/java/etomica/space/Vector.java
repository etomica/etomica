/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;


import com.fasterxml.jackson.databind.annotation.JsonSerialize;
import etomica.math.function.IFunction;
import etomica.meta.annotations.IgnoreProperty;
import etomica.meta.serializers.VectorSerializer;
import etomica.space1d.Vector1D;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;
import etomica.spaceNd.VectorND;
import etomica.util.random.IRandom;

/**
 * Interface containing vector operations, accessor, and mutator methods.
 */
@JsonSerialize(using = VectorSerializer.class)
public interface Vector {

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

    default double[] toArray() {
        double[] arr = new double[this.getD()];
        this.assignTo(arr);
        return arr;
    }

    /**
     * Returns true if all corresponding elements of this and the given vector
     * are equal; returns false otherwise.
     */
    public boolean equals(Vector v);

    /**
     * Vector components corresponding to the given index. For example,
     * x-component is given for i = 0. If index exceeds dimension of vector, no
     * error is given and the last element of the vector is returned.
     */
    public double getX(int i);

    /**
     * Returns the square magnitude of this vector, e.g., x^2 + y^2 for D = 2.
     */
    public double squared();

    /**
     * Returns true if all components of this vector are zero; false otherwise.
     */
    @IgnoreProperty
    public boolean isZero();

    /**
     * Returns the dot product of this vector with the given one.
     */
    public double dot(Vector u);

    /**
     * Returns true if any element of the vector is not-a-number.
     */
    @IgnoreProperty
    public boolean isNaN();

    /**
     * Returns square of vector resulting from subtracting given vector
     * from this one.  Neither vector is changed by this operation.
     */
    public double Mv1Squared(Vector v1);

    /**
     * Sets the i-th component of the vector to the given value d. If index
     * exceeds dimension of vector, no error is given and the last element of
     * the vector is assigned.
     */
    void setX(int i, double d);

    /**
     * Sets the components of this vector equal to those of the given vector.
     */
    void E(Vector u);

    /**
     * Sets all components of this vector equal to the given value.
     */
    void E(double a);

    /**
     * Sets the components of this vector equal to the corresponding elements of
     * the given array. Inverse of the assignTo method.
     */
    void E(double[] a);

    /**
     * Plus-equals (+=) operation. Increments each component of this vector by
     * the corresponding value in the given vector.
     */
    void PE(Vector u);

    /**
     * Plus-equals (+=) operation, causing a constant value to be added to all
     * components of this vector.
     */
    void PE(double a);

    /**
     * Minus-equals (-=) operation. Decrements each component of this vector by
     * the corresponding value in the given vector.
     */
    void ME(Vector u);

    /**
     * Times-equals (*=) operation. Multiplies each component of this vector by
     * the corresponding value in the given vector.
     */
    void TE(Vector u);

    /**
     * Divide-equals (/=) operation. Divides each component of this vector by
     * the corresponding value in the given vector.
     */
    void DE(Vector u);

    /**
     * Multiplies all components by a constant.
     */
    void TE(double a);

    /**
     * Operation (= a1 * v1); sets the components of this to those of the given
     * vector multiplied by the given constant.
     */
    void Ea1Tv1(double a, Vector v1);

    /**
     * Increments (+=) components by (a1 * v1)
     */
    void PEa1Tv1(double a, Vector v1);

    /**
     * Sets the components of this vector equal to (v1 + v2)
     */
    void Ev1Pv2(Vector v1, Vector v2);

    /**
     * Sets the components of this vector equal to (v1 - v2)
     */
    void Ev1Mv2(Vector v1, Vector v2);

    /**
     * Replaces each component of this vector with its value modulo
     * the corresponding component in u
     */
    void mod(Vector u);

    /**
     * Normalizes this vector, so this.squared() == 1.  Divides all
     * components by Math.sqrt(this.squared()).
     */
    void normalize();

    /**
     * Applies the given function to each element of the vector.
     */
    void map(IFunction f);

    /**
     * Sets this vector equal to the cross product of this vector with the
     * given vector (only works for 3D vectors).
     */
    void XE(Vector u);

    /**
     * Assigns this vector to equal a point chosen randomly on the
     * surface of a unit sphere.
     * @param random the random number generator used to perform the operation
     */
    void setRandomSphere(IRandom random);

    /**
     * Assigns each component to (its own) random value between -0.5 and + 0.5.
     * @param random the random number generator used to perform the operation
     */
    void setRandomCube(IRandom random);

    /**
     * Assigns this vector to equal a point chosen randomly in the volume
     * of a unit sphere.
     * @param random the random number generator used to perform the operation
     */
    void setRandomInSphere(IRandom random);

    void nearestImage(Vector dimensions);


    /**
     * Returns a Vector initialized to the given set of values.
     * Spatial dimension of the Vector is determined by the number of components.
     *
     * @throws IllegalArgumentException if number of components is not greater than 0
     */
    static Vector of(double... components) {
        switch (components.length) {
            case 0:
                throw new IllegalArgumentException("Number of vector components must be greater than 0");
            case 1:
                return new Vector1D(components);
            case 2:
                return new Vector2D(components);
            case 3:
                return new Vector3D(components);
            default:
                return new VectorND(components);

        }
    }

    /**
     * Returns a Vector initialized to the given set of values (cast to double).
     * Spatial dimension of the Vector is determined by the number of components.
     */
    static Vector of(int... components) {
        double[] a = new double[components.length];
        for (int i = 0; i < components.length; i++) {
            a[i] = (double) components[i];
        }
        return Vector.of(a);
    }

    /**
     * Creates a Vector of the given dimension. Uses specific vector types for dimensions 1, 2, and 3, and VectorND for
     * any higher dimension.
     *
     * @param dimension greater than 0
     * @return a Vector of the given dimension.
     */
    static Vector d(int dimension) {
        if (dimension <= 0) {
            throw new IllegalArgumentException("Vector dimension must be greater than 0");
        }
        switch (dimension) {
            case 1:
                return new Vector1D();
            case 2:
                return new Vector2D();
            case 3:
                return new Vector3D();
            default:
                return new VectorND(dimension);
        }

    }

    Vector duplicate();
}
