/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;


/**
 * Basic vector interface containing methods needed for basic vector operations
 * and other accessor methods.
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
    public double getX(int i);

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
}