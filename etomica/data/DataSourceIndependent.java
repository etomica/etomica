package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;

/**
 * Interface for an object that can provide X Data objects on request.  This
 * provides the independent data part of a DataFunction, which downstream data
 * elements can retrieve when they need it.
 */ 
public interface DataSourceIndependent {

    /**
     * Returns the X data for the given dimension
     */
	public DataDoubleArray getIndependentData(int i);

    /**
     * Returns the DataInfo for the given dimension
     */
    public DataInfoDoubleArray getIndependentDataInfo(int i);
    
    /**
     * Returns the number of independent data dimensions
     */
    public int getIndependentArrayDimension();
}