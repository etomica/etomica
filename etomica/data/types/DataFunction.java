package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;
import etomica.utility.Function;


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
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

/*
 * History
 * Created on Jun 16, 2005 by kofke
 */
public class DataFunction extends Data implements DataArithmetic {

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
    public DataFunction(DataDoubleArray[] independentData, DataDoubleArray dependentData) {
        super(new DataInfo(dependentData.getDataInfo().getLabel(), 
                dependentData.getDataInfo().getDimension(), getFactory(independentData)));
        if(dependentData.getArrayDimension() != independentData.length) {
            throw new IllegalArgumentException("Dimension of dependent data is not compatible with number of independent data elements");
        }
        for(int i=0; i<independentData.length; i++) {
            if(independentData[i].getArrayDimension() != 1) {
                throw new IllegalArgumentException("All independent data must be of dimension 1");
            }
        }
        this.dependentData = dependentData;
        this.independentData = (DataDoubleArray[])independentData.clone();
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
        this(independentData, makeDependentData(label, dimension, independentData));
    }
    
    //used by constructor above
    private static DataDoubleArray makeDependentData(String label, Dimension dimension, 
            DataDoubleArray[] independentData) {
        int[] size = new int[independentData.length];
        for (int i=0; i<size.length; i++) {
            size[i] = independentData[i].getLength();
        }
        return new DataDoubleArray(label,dimension,size);
    }

    /**
     * Copy constructor. New DataFunction instance and original share the same 
     * independent data instances (and DataInfo); new and original have different
     * instances of dependent data.
     */
    public DataFunction(DataFunction data) {
        super(data);
        this.independentData = data.independentData;
        dependentData = new DataDoubleArray(data.dependentData);
    }
    
    /**
     * Returns a copy of this instance, holding the same instances of the independent
     * data and a new instance of the dependent data.
     */
    public Data makeCopy() {
        return new DataFunction(this);
    }


    /**
     * Applies the E (equals) method to the dependent data. 
     * 
     * @throws ClassCastException if the given data is not an instance of DataFunction.
     * 
     */
    public void E(Data data) {
        dependentData.E(((DataFunction)data).dependentData);
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
     * Returns the dependent data.
     */
    public DataDoubleArray getYData() {
        return dependentData;
    }
    
    /**
     * Returns the dimension of the independent data, the number of
     * values the function depends upon.
     */
    public int getXDimension() {
        return independentData.length;
    }
    
    /**
     * Plus-equals (+=) operation. Element-by-element addition of the values in
     * the given array to those in this one.
     * 
     * @throws ArrayIndexOutOfBoundException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void PE(DataArithmetic y) {
        if (y instanceof DataFunction) {
            y = ((DataFunction)y).dependentData;
        }
        dependentData.PE(y);
    }

    /**
     * Minus-equals (-=) operation. Element-by-element subtraction of the values
     * in the given array from those in this one.
     * 
     * @throws ArrayIndexOutOfBoundException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void ME(DataArithmetic y) {
        if (y instanceof DataFunction) {
            y = ((DataFunction)y).dependentData;
        }
        dependentData.ME(y);
    }

    /**
     * Times-equals (*=) operation. Applied element-by-element.
     * 
     * @throws ArrayIndexOutOfBoundException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void TE(DataArithmetic y) {
        if (y instanceof DataFunction) {
            y = ((DataFunction)y).dependentData;
        }
        dependentData.TE(y);
    }

    /**
     * Divide-equals (/=) operation. Applied element-by-element.
     * 
     * @throws ArrayIndexOutOfBoundException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void DE(DataArithmetic y) {
        if (y instanceof DataFunction) {
            y = ((DataFunction)y).dependentData;
        }
        dependentData.DE(y);
    }

    /**
     * Equals (=) operation. Sets all elements to given value.
     */
    public void E(double y) {
        dependentData.E(y);
    }

    /**
     * Plus-equals (+=) operation. Increments all elements by the given value.
     */
    public void PE(double y) {
        dependentData.PE(y);
    }

    /**
     * Times-equals (*=) operation. Multiplies all elements by the given value.
     */
    public void TE(double y) {
        dependentData.TE(y);
    }

    /**
     * Replaces all values by the value of the function applied to each.
     */
    public void map(Function function) {
        dependentData.map(function);
    }

    /**
     * Returns the length of the wrapped array.
     */
    public int getLength() {
        return dependentData.getLength();
    }
    
    /**
     * Returns true if any of the dependent data are NaN.
     */
    public boolean isNaN() {
        return dependentData.isNaN();
    }
    
    /**
     * Assigns dependent-data values to the given array.
     */
    public void assignTo(double[] array) {
        dependentData.assignTo(array);
    }
    
    /**
     * Returns the i-th dependentData value.
     */
    public double getValue(int i) {
        return dependentData.getValue(i);
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
        string.append(dependentData.toString());
        return string.toString();
    }
    
    private final DataDoubleArray[] independentData;
    private final DataDoubleArray dependentData;
    
    /**
     * Returns a factory for a DataFunction. Each new instance made by the factory
     * will be a DataFunction having its own copy of the given dependent
     * data, while sharing the same instances of the given independent data.
     */
    public static DataFactory getFactory(DataDoubleArray[] independentData) {
        return new Factory(independentData);
    }

    /**
     * DataFactory that constructs DataFunction instances all having the same
     * set of independent data.  Instantiate using the static DataFunction getFactory method.
     */
    public static class Factory implements DataFactory, Serializable {
        
        final DataDoubleArray[] independentData;
        
        Factory(DataDoubleArray[] independentData) {
            this.independentData = (DataDoubleArray[])independentData.clone();
        }
        
        /**
         * Makes a new DataFunction with the given label and dimension for the
         * dependent data, and with the prototype instances of the independent data.
         */
        public Data makeData(String label, Dimension dimension) {
            return new DataFunction(label, dimension, independentData);
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
        
        /**
         * Returns an array of the lengths of the independent-data arrays
         * in the DataFunction made by this factory.
         */
        public int[] getIndependentDataSizes() {
            int[] size = new int[independentData.length];
            for(int i=0; i<size.length; i++) {
                size[i] = independentData[i].getLength();
            }
            return size;
        }
        
    }//end of Factory

}
