package etomica.data.types;

import etomica.Data;
import etomica.utility.Function;

/**
 * Interface for a Data class that supports arithmetic operations on its data.
 * Methods throw IllegalArgumentException if this instance and operand are not
 * compatible.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jun 15, 2005 by kofke
 */
public interface DataArithmetic extends DataNumeric {

    /**
     * Equals (=) operation. Sets the data in this instance equal to the values
     * in the given object.
     */
    public void E(Data y);

    /**
     * Plus-equals (+=) operation.
     */
    public void PE(DataArithmetic y);

    /**
     * Minus-equals (-=) operation.
     */
    public void ME(DataArithmetic y);

    /**
     * Times-equals (*=) operation.
     */
    public void TE(DataArithmetic y);

    /**
     * Divide-equals (/=) operation.
     */
    public void DE(DataArithmetic y);

    /**
     * Equals (=) operation, sets all values in data equal to the given value.
     */
    public void E(double y);

    /**
     * Plus-equals (+=) operation, adding given value to all values in data.
     * 
     * @param y
     */
    public void PE(double y);

    /**
     * Times-equals (*=) operation, multiplying all values in data by given
     * value.
     * 
     * @param y
     */
    public void TE(double y);

    /**
     * Maps the function on all data values, replace each with the value given
     * by the function applied to it.
     */
    public void map(Function function);

    /**
     * Returns the number of values held by the data instance.
     */
    public int getLength();//TODO consider changing to setNValues

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
     * Returns a new array formed from the values held by the data instance.
     */
    public double[] toArray();

    /**
     * Returns true if any data value is true for Double.isNaN
     */
    public boolean isNaN();
}