/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.veos;

public interface IBigValue {
    
    /**
     * Returns an instance of the class representing the value zero.
     */
    public IBigValue getZero();

    /**
     * Returns true if this instance has a value of zero.
     */
    public boolean isZero();

    /**
     * Returns the value held by the instance (in true form, not as log, and with sign applied).
     */
    public double value();
    
    
    /**
     * Returns true if this instance has a value greater than zero.
     */
    public boolean isPositive();

    /**
     * Returns logarithm of the value represented by the instance. Returns NaN if the
     * value is negative.
     */
    public double lnValue();

    /**
     * Equals (=) operation.  Assigns the given value to this instance, replacing current value.
     */
    public void E(double value);

    /**
     * Equals (=) operation. Assigns given value (specified via another instance of LogValue) to this instance, 
     * replacing current value.
     */
    public void E(IBigValue v);

    /**
     * Plus-equals (+=) operation.  Replaces the current value with the one obtained by adding
     * the given value, considering signs of both.
     */
    public void PE(IBigValue x);

    /**
     * Plus-equals (+=) a1 times v1.
     */
    public void PEa1Tv1(double a1, IBigValue v1);

    /**
     * Times-equals (*=) operation.  Replaces the current value with the one obtained by multiplying
     * the given value, considering signs of both.
     */
    public void TE(IBigValue x);

    /**
     * Divide-equals (/=) operation.  Replaces the current value with the one obtained by dividing
     * the given value, considering signs of both.
     */
    public void DE(IBigValue x);

    /**
     * Times-equals (*=). Replaces current value by multiplying it by the given value.
     * @param a1
     */
    public void TE(double a1);

    /**
     * Times-equals (*=) operation, multiplying by a value given as its logarithm.  Value (i.e., the value
     * whose log is given as the argument) is assumed to be positive.
     */
//    public void TEln(double lnx);

}