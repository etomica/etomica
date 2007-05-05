package etomica.util;

import etomica.data.Data;
import etomica.data.DataInfo;

/**
 * Interface for the basic features of a general function that maps an object
 * into Data.
 */
public interface FunctionGeneral {

    /**
     * Method that performs the mapping defined by the function.
     * 
     * @param x
     *            the input object, typically a Double, Vector, Tensor, Atom,
     *            etc.
     * @return the output instance of Data defined by this function
     */
    public Data f(Object x);

    /**
     * Returns a DataInfo that describes the Data object that is returned by this function's f method.
     */
    public DataInfo getDataInfo();

}
