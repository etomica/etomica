package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;

/**
 * Data object that wraps two arrays of double, such that one is considered to
 * have a dependence on the other. The dependence is not enforced by this data
 * structure, rather this structure is used to encapsulate a set of x-t value
 * pairs, where the function x(t) is applied by some other means.  Data arrays
 * are accessed via the getData method (for x) and getTData (for t, which returns
 * a DataDoubleArray instance).
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

/*
 * History Created on Jun 15, 2005 by kofke
 */
public class DataFunction extends DataDoubleArray {

    /**
     * Function is x(t).
     * 
     * @param xDataInfo
     *            dataInfo for the dependent variable
     * @param tDataInfo
     *            dataInfo for the independent variable
     */
    public DataFunction(DataInfo xDataInfo, DataInfo tDataInfo) {
        super(xDataInfo);
        t = new DataDoubleArray(tDataInfo);
    }

    /**
     * Copy constructor.
     */
    public DataFunction(DataFunction data) {
        super(data);
        t = (DataDoubleArray) data.t.makeCopy();
    }

    /**
     * Returns a copy of this instance. Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataFunction(this);
    }

    /**
     * Copies both arrays in the given object to the arrays in this one.
     * 
     * @throws ClassCastException
     *             if the given Data object is not of type DataFunction
     */
    public void E(Data y) {
        super.E(y);
        t.E(((DataFunction) y).t);
    }

    /**
     * Sets the length of both wrapped arrays to the given value, erasing
     * previous values, if the given length is different from the current
     * length.
     * 
     * @throws IllegalArgumentException
     *             if n < 0
     */
    public void setLength(int n) {
        super.setLength(n);
        t.setLength(n);
    }

    /**
     * Returns the data object wrapping the independent-variable array.
     */
    public DataDoubleArray getTData() {
        return t;
    }

    /**
     * Returns a string formed from the dataInfo label and the two wrapped
     * arrays.
     */
    public String toString() {
        return dataInfo.getLabel() + " (" + super.toString() + "), ("
                + t.toString() + ")";
    }

    private final DataDoubleArray t;
}