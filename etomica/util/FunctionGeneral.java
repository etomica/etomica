package etomica.util;

import etomica.data.IData;
import etomica.data.IDataInfo;

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
    public IData f(Object x);

    /**
     * Returns a DataInfo that describes the Data object that is returned by this function's f method.
     */
    public IDataInfo getDataInfo();

}
