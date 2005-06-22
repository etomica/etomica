package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;

/**
 * Data object that wraps two arrays of double, such that one is considered
 * to have a dependence on the other.  The dependence is not enforced
 * by this data structure, rather this structure is used to encapsulate a
 * set of x-t value pairs, where the function x(t) is applied by some other
 * means.
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
     * @param xDataInfo  dataInfo for the dependent variable
     * @param tDataInfo  dataInfo for the independent variable
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
        t = (DataDoubleArray)data.t.makeCopy();
    }
    
    /**
     * Returns a copy of this instance.  Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataFunction(this);
    }

    public void E(Data y) {
        super.E(y);
        t.E(((DataFunction)y).t);
    }

    public void setLength(int n) {
        super.setLength(n);
        t.setLength(n);
    }

    public DataDoubleArray getTData() {
        return t;
    }
    
    public String toString() {
        return dataInfo.getLabel() + " ("+super.toString()+"), ("+t.toString()+")";
    }

    private final DataDoubleArray t;
}