/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.math.function.IFunction;

/**
 * Abstract container of simulation data and information describing it. Data
 * objects may encapsulate any type of data (such as one or more double values),
 * and are transmitted from a DataSource to a DataSink, perhaps passing through
 * a series of data processing elements along the way. The abstract Data class
 * provides methods to obtain descriptive information about the data (its
 * DataInfo), a method to copy the data to another existing data object
 * instance, and a method to make a copy of the data to a new data object
 * instance. Different subclasses of Data are generally not interchangable, and
 * objects processing or receiving Data instances often expect a specific
 * subclass.
 * <p>
 * The data held by a Data instance is mutable, but the structure of a Data
 * instance is not. The "structure" of a Data instance refers to the number and
 * type of data it holds. A DataDoubleArray, for example, holds an array of
 * double; the values in the array may be changed, but the length of the array
 * (the structure) cannot. Thus any object holding a reference to a Data
 * instance can be assured that the number and type of data it holds will not be
 * changed by something else that also references the instance.
 * <p>
 * Information about the data is encapsulated in a DataInfo instance that is
 * held by the Data instance. The DataInfo holds a descriptive label, a
 * Dimension instance that describes the physical dimensions (e.g., mass, time)
 * of the data, and a DataFactory that will make new independent instances of
 * the same Data type and with the same structure. The DataInfo instance is,
 * like the Data structure, completely immutable.
 * 
 * @author David Kofke and Andrew Schultz
 * 
 * @see DataInfo
 */

public interface IData {

    /**
     * Deep-copies the data from the given object to this one.
     * 
     * @throws ClassCastException
     *             if the given data object is not of a type that can be copied
     *             to this one.
     */
    public void E(IData data);
    
    /**
     * Plus-equals (+=) operation.
     */
    public void PE(IData y);

    /**
     * Minus-equals (-=) operation.
     */
    public void ME(IData y);

    /**
     * Times-equals (*=) operation.
     */
    public void TE(IData y);

    /**
     * Divide-equals (/=) operation.
     */
    public void DE(IData y);

    /**
     * Equals (=) operation, sets all values in data equal to the given value.
     */
    public void E(double y);

    /**
     * Plus-equals (+=) operation, adding given value to all values in data.
     */
    public void PE(double y);

    /**
     * Times-equals (*=) operation, multiplying all values in data by given
     * value.
     */
    public void TE(double y);

    /**
     * Divide-equals (*=) operation, dividing all values in data by given
     * value.
     */
    public void DE(double y);

    /**
     * Maps the function on all data values, replace each with the value given
     * by the function applied to it.
     */
    public void map(IFunction function);

    /**
     * Returns the number of values held by the data instance.
     */
    public int getLength();//TODO consider changing to getNValues

    /**
     * Returns the i-th data value.
     * 
     * @param i
     *            index of the desired data value
     * @return the data value
     * @throws IllegalArgumentException
     *             if i >= getLength()
     */
    public double getValue(int i);

    /**
     * Fills the given array with the values held by this Data instance.
     * 
     * @throws IllegalArgumentException
     *             if the length of the array is not equal to getLength
     */
    public void assignTo(double[] array);

    /**
     * Returns true if any data value is true for Double.isNaN
     */
    public boolean isNaN();
    
}
