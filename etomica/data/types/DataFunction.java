package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.units.Dimension;


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
     * @param yData
     *            the data defined on the space of independent data; a
     *            double[] with getDimension equal to
     *            independentData.length
     * 
     * @throws IllegalArgumentException
     *             if length of independentData array does not equal
     *             dependentData.getArrayDimension(), or if any of the
     *             independent data have array dimension not equal to 1
     */
    public DataFunction(DataDoubleArray[] independentData, double[] yData) {
        super(getArrayShape(independentData), yData);
        this.independentData = (DataDoubleArray[])independentData.clone();
    }
    
    /**
     * Forms a DataFunction using the given independent data instances.  The 
     * dependent data is constructed using the array dimensions from the 
     * indepdendent data.
     * 
     * @param independentData
     *            array of DataDoubleArray, each representing independent-data
     *            values in one dimension. Domain of data is determined by this
     *            array of independent data. Given array is cloned for internal
     *            use, while instances in array are used directly.
     * 
     * @throws IllegalArgumentException
     *            if any of the independent data have array dimension not 
     *            equal to 1
     */
    public DataFunction(DataDoubleArray[] independentData) {
        super(getArrayShape(independentData));
        this.independentData = (DataDoubleArray[])independentData.clone();
    }
    
    protected static int[] getArrayShape(DataDoubleArray[] independentArray) {
        int[] arrayShape = new int[independentArray.length];
        for (int i=0; i<arrayShape.length; i++) {
            if (independentArray[i].getArrayDimension() != 1) {
                throw new IllegalArgumentException("Array dimension must be 1");
            }
            arrayShape[i] = independentArray[i].getLength();
        }
        return arrayShape;
    }
    
    /**
     * Forms a DataFunction with the given array shape.
     */
    public DataFunction(int[] arrayShape) {
        super(arrayShape);
        independentData = new DataDoubleArray[arrayShape.length];
        for(int i=0; i<arrayShape.length; i++) {
            this.independentData[i] = new DataDoubleArray(arrayShape[i]);
        }
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
        StringBuffer string = new StringBuffer();
        for(int i=0; i<independentData.length; i++) {
            string.append("\n"); //newline?
            string.append(independentData.toString());
        }
        string.append(toString());
        return string.toString();
    }
    
    private final DataDoubleArray[] independentData;//shadows the field in Factory
    
    public static class DataInfoFunction extends DataInfoDoubleArray {
        public DataInfoFunction(String label, Dimension dimension, DataInfoDoubleArray[] independentInfo) {
            super(label, dimension, getArrayShape(independentInfo));
            this.independentInfo = (DataInfo[])independentInfo.clone();
        }
        
        private static int[] getArrayShape(DataInfoDoubleArray[] independentInfo) {
            int[] arrayShape = new int[independentInfo.length];
            for (int i=0; i<arrayShape.length; i++) {
                arrayShape[i] = independentInfo[i].getArrayShape()[0];
            }
            return arrayShape;
        }
        
        public DataInfo getIndependentInfo(int i) {
            return independentInfo[i];
        }
        
        protected final DataInfo[] independentInfo;
    }
    
    public static class DataInfoFunctionFactory extends DataInfoDoubleArrayFactory {
        protected DataInfoFunctionFactory(DataInfoFunction template) {
            super(template);
            independentInfo = (DataInfoFunction[])template.independentInfo.clone();
        }
        
        public DataInfo makeDataInfo() {
            return new DataInfoFunction(label, dimension, independentInfo);
        }
        
        /**
         * Sets array of independent DataInfo.  The array is copied so further
         * changes made to the given array will not be affect this factory.
         */
        public void setIndependentInfo(DataInfoDoubleArray[] newIndependentInfo) {
            independentInfo = (DataInfoFunction[])newIndependentInfo.clone();
        }
        
        /**
         * Returns a copy of array of independent DataInfo.
         */
        public DataInfoDoubleArray[] getIndependentInfo() {
            return (DataInfoDoubleArray[])independentInfo.clone();
        }
        
        protected DataInfoDoubleArray[] independentInfo;
    }
}
