package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.units.Dimension;
import etomica.units.Null;


/**
 * Collects two or more DataDoubleArray instances, and organizes them into
 * "dependent" and "independent" subgroups. Normally these data represent a
 * functional dependence, in which the dependent data were calculated at each
 * point in the domain of the independent data. However, nothing about this
 * class enforces the dependence, it merely classifies the data into the two
 * groups.
 * <p>
 * Multidimensional functions hold the independent data as an array of
 * one-dimensional DataDoubleArrays, and the dependent data is represented by a
 * single DataDoubleArray that is formed with a corresponding number of
 * dimensions.
 * <p>
 * The DataInfo instance for the DataFunction has a the same label and dimension
 * as the that for the dependent data (although their factories differ).
 * <p>
 * All arithmetic operations apply only to the dependent data. Independent data
 * values are unaffected by them.
 * <p>
 * Note that all instances created by the same factory will hold the same instances
 * of the independent data.  Thus if the independent data is changed in any instance,
 * it will be reflected in all instances that were constructed by the factory.
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

/*
 * History
 * Created on Jun 16, 2005 by kofke
 */
public class DataFunction extends DataDoubleArray {

    /**
     * Forms a DataFunction using the given independent and dependent data
     * instances.
     * 
     * @param independentData
     *            array of DataDoubleArray, each representing independent-data
     *            values in one dimension. Domain of data is determined by this
     *            array of independent data. Given array is cloned for internal
     *            use, while instances in array are used directly.
     * @param dependentData
     *            the data defined on the space of independent data; a
     *            DataDoubleArray with getDimension equal to
     *            independentData.length
     * 
     * @throws IllegalArgumentException
     *             if length of independentData array does not equal
     *             dependentData.getArrayDimension(), or if any of the
     *             independent data have array dimension not equal to 1
     */
    public DataFunction(String label, Dimension dimension, DataDoubleArray[] independentData, double[] yData) {
        super(label, dimension, getArrayShape(independentData), yData, getFactory(independentData));
        this.independentData = (DataDoubleArray[])independentData.clone();
    }
    
    protected static int[] getArrayShape(DataDoubleArray[] independentData) {
        int[] arrayShape = new int[independentData.length];
        for (int i=0; i<arrayShape.length; i++) {
            arrayShape[i] = independentData[i].getLength();
        }
        return arrayShape;
    }
    
    /**
     * Forms a DataFunction with the given label and dimension, and with a
     * domain of data defined by the array of independent data. The array of
     * dependent data is constructed with shape given by the lengths of the
     * independent data array.
     * 
     * @throws IllegalArgumentException
     *             if any of the independent data have array dimension not equal
     *             to 1
     */
    public DataFunction(String label, Dimension dimension, DataDoubleArray[] independentData) {
        super(label, dimension, (DataFunction.Factory)getFactory(independentData));
        this.independentData = (DataDoubleArray[])independentData.clone();
    }
    
    /**
     * Copy constructor. New DataFunction instance and original share the same 
     * independent data instances (and DataInfo); new and original have different
     * instances of dependent data.
     */
    public DataFunction(DataFunction data) {
        super(data);
        this.independentData = data.independentData;
    }
    
    /**
     * Returns a copy of this instance, holding the same instances of the independent
     * data and a new instance of the dependent data.
     */
    public Data makeCopy() {
        return new DataFunction(this);
    }


    /**
     * Returns the i-th set of independent data.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if i < 0 or i >= getXDimension()
     */
    public DataDoubleArray getXData(int i) {
        return independentData[i];
    }

    /**
     * Returns the dimension of the independent data, the number of
     * values the function depends upon.
     */
    public int getXDimension() {
        return independentData.length;
    }
    
    /**
     * Returns a string formed from the encapsulated data objects.
     */
    public String toString() {
        StringBuffer string = new StringBuffer(dataInfo.getLabel());
        for(int i=0; i<independentData.length; i++) {
            string.append("\n"); //newline?
            string.append(independentData.toString());
        }
        string.append(toString());
        return string.toString();
    }
    
    private final DataDoubleArray[] independentData;//shadows the field in Factory
    
    /**
     * Returns a factory for a DataFunction. Each new instance made by the factory
     * will be a DataFunction having its own copy of the given dependent
     * data, while sharing the same instances of the given independent data
     * (not the array itself, but the elements it holds are common to all
     * DataFunction instances made by the factory).
     */
    public static DataFactory getFactory(DataDoubleArray[] independentData) {
        return new Factory(independentData);
    }

    /**
     * DataFactory that constructs DataFunction instances all having the same
     * set of independent data.  Instantiate using the static DataFunction getFactory method.
     */
    public static class Factory extends DataDoubleArray.Factory {
        
        final DataDoubleArray[] independentData;
        
        Factory(DataDoubleArray[] independentData) {
            super(DataFunction.getArrayShape(independentData));
            this.independentData = (DataDoubleArray[])independentData.clone();
        }
        
        /**
         * Makes a new DataFunction with the given label and dimension for the
         * dependent data, and with the prototype instances of the independent data.
         */
        public Data makeData(String label, Dimension dimension) {
            return new DataFunction(label, dimension, (DataDoubleArray[])independentData.clone());
        }
        
        /**
         * Returns DataFunction.class, indicating that this factory makes DataFunction instances.
         */
        public Class getDataClass() {
            return DataFunction.class;
        }
        
        /**
         * Returns a clone of the independent-data array.
         */
        public DataDoubleArray[] getIndependentData() {
            return (DataDoubleArray[])independentData.clone();
        }
        
    }//end of Factory

}
