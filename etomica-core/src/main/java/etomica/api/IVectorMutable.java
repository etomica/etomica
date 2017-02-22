/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;


/**
 * Vector interface containing methods needed for vector operations which alter
 * the contents of the vector.
 */
public interface IVectorMutable extends IVector {

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
    public void PE(double a);

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
     * Normalizes this vector, so this.squared() == 1.  Divides all 
     * components by Math.sqrt(this.squared()).
     */
    public void normalize();

    /**
     * Applies the given function to each element of the vector.
     */
    public void map(IFunction f);
    

    /**
     * Sets this vector equal to the cross product of this vector with the
     * given vector (only works for 3D vectors).
     */
    public void XE(IVector u);
}