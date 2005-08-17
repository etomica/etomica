package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;
import etomica.utility.Function;


/**
 * Collects two or more Data instances, and organizes them into "dependent"
 * and "independent" subgroups. Normally these data represent a functional 
 * dependence, in which the dependent data were calculated at each point in the
 * domain of the independent data. However, nothing about this class enforces
 * the dependence, it merely classifies the data into the two groups.
 * <p>
 * Data sinks typically need to manipulate DataFunction instances to extract the individual
 * data elements in them.
 * <p>
 * The DataInfo instance for the DataFunction has a the same label and dimension as the 
 * that for the dependent data (although their factories differ).
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
     *            array of independent data.
     * @param dependentData
     *            the data defined on the space of independent data; a
     *            DataDoubleArray with getDimension equal to
     *            independentData.length
     */
    public DataFunction(DataDoubleArray[] independentData, DataDoubleArray dependentData) {
        super(new DataInfo(dependentData.getDataInfo().getLabel(), 
                dependentData.getDataInfo().getDimension(), getFactory(independentData)));
        this.dependentData = dependentData;
        this.independentData = (DataDoubleArray[])independentData.clone();
    }
    
    /**
     * Forms a DataFunction with the given label and dimension, and with a
     * domain of data defined by the array of independent data.  The array
     * of dependent data is constructed with shape given by the lengths of
     * the independent data array.
     */
    public DataFunction(String label, Dimension dimension, 
            DataDoubleArray[] independentData) {
        super(new DataInfo(label, dimension, getFactory(independentData)));
        int[] size = new int[independentData.length];
        for (int i=0; i<size.length; i++) {
            size[i] = independentData[i].getLength();
        }
        dependentData = new DataDoubleArray(label,dimension,size);
        this.independentData = (DataDoubleArray[])independentData.clone();
    }

    /**
     * Copy constructor.  Performs a deep copy, forming a new
     * DataGroup with new instances of copies of this instance's
     * data objects.
     */
    public DataFunction(DataFunction data) {
        super(data);
        independentData = new DataDoubleArray[data.getXDimension()];
        for (int i=0; i<independentData.length; i++) {
            independentData[i] = new DataDoubleArray(data.getXData(i));
        }
        dependentData = new DataDoubleArray(data.getYData());
    }
    
    /**
     * Returns a deep copy of this instance.  Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.  All data
     * objects in this group are copied.
     */
    public Data makeCopy() {
        return new DataFunction(this);
    }


    /**
     * Applies the E method to all data elements held, in a one-to-one
     * correspondence with the elements in the given data group. 
     * 
     * @throws ClassCastException if the given data is not an instance of DataFunction.
     * @throws IllegalArgumentException if the given DataGroup has a different number of data elements than this DataGroup.
     * 
     */
    public void E(Data data) {
        dependentData.E(((DataFunction)data).dependentData);
    }
    
    /**
     * Returns the i-th set of independent data.
     * @param i
     * @return
     */
    public DataDoubleArray getXData(int i) {
        return independentData[i];
    }

    public DataDoubleArray getYData() {
        return dependentData;
    }
    
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
    
    public boolean isNaN() {
        return dependentData.isNaN();
    }
    
    public void assignTo(double[] array) {
        dependentData.assignTo(array);
    }
    
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
    
    final DataDoubleArray[] independentData;
    final DataDoubleArray dependentData;
    
    /**
     * Returns a factory for a data group holding copies of the
     * given data objects.  Each new instance made by the factory
     * will be a DataGroup having its owns copy of the given data.
     * When the factory is made, the given Data array is cloned, 
     * but the data instances it holds are not.
     */
    public static DataFactory getFactory(DataDoubleArray[] data) {
        return new Factory(data);
    }

    public static class Factory implements DataFactory, Serializable {
        
        final DataDoubleArray[] independentData;
        
        Factory(DataDoubleArray[] independentData) {
            this.independentData = independentData;
        }
        
        public Data makeData(String label, Dimension dimension) {
            return new DataFunction(label, dimension, independentData);
        }
        
        public Class getDataClass() {
            return DataFunction.class;
        }
        
        /**
         * Returns a clone of the independentData array.
         */
        public DataDoubleArray[] getIndependentData() {
            return (DataDoubleArray[])independentData.clone();
        }
        
        public int[] getIndependentDataSizes() {
            int[] size = new int[independentData.length];
            for(int i=0; i<size.length; i++) {
                size[i] = independentData[i].getLength();
            }
            return size;
        }
        
    }

}
