package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;


/**
 * Gathers one or more Data instances into a single Data object.  The primary
 * use is to enable a data handler to take a single data stream and emit multiple
 * different streams from it.  An example is given by the accumulators, which
 * take data and compile one or more new sets of data from it, and send them out
 * bundled in a DataGroup.
 * <p>
 * The Data instances held by the group are set at construction and cannot be changed.
 * <p>
 * Data sinks typically need to manipulate DataGroups to extract the individual
 * data elements in them.
 * <p>
 * The DataInfo instance for the DataGroup has a label that describes the collection,
 * and a dimension that is the dimension of all of the data objects it holds (if they
 * are all the same) or MIXED (if they are not).
 *
 * @author David Kofke and Andrew Schultz
 *
 */

/*
 * History
 * Created on Jun 16, 2005 by kofke
 */
public class DataGroup extends Data {

    /**
     * Forms a data group from the given array of data objects. Given data array
     * is cloned, but the data instances it holds are kept by this DataGroup.
     * <p>
     * Dimension associated with the DataGroup is Dimension.MIXED if the Data it
     * holds have different dimensions; otherwise dimension is the common
     * dimension of the grouped Data instances.
     * 
     * @param label
     *            a descriptive label for the data collection
     * @param data
     *            array of data to be encapsulated in this group
     * 
     * @throws NullPointerException
     *             if any of the elements of the Data array are null
     */
    public DataGroup(Data[] data) {
        super();
        this.data = (Data[])data.clone();
    }
    
    /**
     * Copy constructor.  Performs a deep copy, forming a new
     * DataGroup with new instances of copies of this instance's
     * data objects.
     */
    public DataGroup(DataGroup group) {
        super();
        this.data = new Data[group.data.length];
        for(int i=0; i<data.length; i++) {
            if(group.data[i] != null) {
                this.data[i] = group.data[i].makeCopy();
            }
        }
    }
    
    /**
     * Returns a deep copy of this instance.  Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.  All data
     * objects in this group are copied.
     */
    public Data makeCopy() {
        return new DataGroup(this);
    }


    /**
     * Applies the E method to all Data elements held, in a one-to-one
     * correspondence with the elements in the given data group.
     * 
     * @throws ClassCastException
     *             if the given data is not an instance of DataGroup.
     * @throws IllegalArgumentException
     *             if the given DataGroup has a different number of data
     *             elements than this DataGroup.
     *  
     */
    public void E(Data newData) {
        if(((DataGroup)newData).data.length != data.length) {
            throw new IllegalArgumentException("Attempt to copy data groups of different length: (this.length, argument's length): ("+this.data.length+", "+((DataGroup)newData).data.length+")");
        }
        for (int i=0; i<this.data.length; i++) {
            data[i].E(((DataGroup)newData).getData(i));
        }
    }

    /**
     * Returns the i-th data element in the group, counting from 0.
     * 
     * @throws ArrayIndexOutOfBounds
     *             exception if the given value does not reference a legitimate
     *             element
     */
    public Data getData(int i) {
        return data[i];
    }
    
    /**
     * Returns the number of data objects encapsulated by this DataGroup.
     */
    public int getNData() {
        return data.length;
    }

    /**
     * Returns a string formed from the encapsulated data objects.
     */
    public String toString() {
        StringBuffer string = new StringBuffer("");
        for(int i=0; i<data.length; i++) {
            string.append("\n"); //newline?
            string.append(data.toString());
        }
        return string.toString();
    }
    
    protected final Data[] data;
    
    /**
     * Returns a factory for a data group holding copies of the
     * given DataFactory objects.  Each new instance made by the factory
     * will be a DataGroup having data made by the given factories.
     * When the factory is made, the given DataFactory array is cloned, 
     * but the DataFactory instances it holds are not.
     */
    public static Factory getFactory(DataFactory[] dataFactory) {
        return new Factory(dataFactory);
    }

    /**
     * DataFactory that makes DataGroups with copies of a specific set of Data.
     * Instantiate using the static DataGroup.getFactory method.
     */
    public static class Factory implements DataFactory, Serializable {
        
        final DataFactory[] dataFactory;
        
        Factory(DataFactory[] dataFactory) {
            this.dataFactory = (DataFactory[])dataFactory.clone();
        }
        
        /**
         * Makes a new DataGroup from the Data in the prototype DataGroup.
         * Dimension is given to adhere to DataFactory interface, but it is not used.
         */
        public Data makeData() {
            Data[] newData = new Data[dataFactory.length];
            for(int i=0; i<dataFactory.length; i++) {
                newData[i] = dataFactory[i].makeData();
            }
            //drop dimension
            return new DataGroup(newData);
        }

        /**
         * Returns a clone of the array of Data held by the prototype DataGroup.
         */
        public DataFactory[] getDataFactory() {
            return (DataFactory[])dataFactory.clone();
        }

        /**
         * Returns DataGroup.class, indicating that this DataFactory produces a DataGroup.
         */
        public Class getDataClass() {
            return DataGroup.class;
        }
    }

    public static class DataInfoGroup extends DataInfo {
        public DataInfoGroup(String label, Dimension dimension, DataInfo[] subDataInfo) {
            super(label, dimension, DataGroup.getFactory(makeFactories(subDataInfo)));
            this.subDataInfo = (DataInfo[])subDataInfo.clone();
        }
        
        protected DataInfoGroup(String label, Dimension dimension, DataInfo[] subDataInfo, DataGroup.Factory factory) {
            super(label, dimension, factory);
            this.subDataInfo = (DataInfo[])subDataInfo.clone();
        }
        
        private static DataFactory[] makeFactories(DataInfo[] subDataInfo) {
            DataFactory[] factories = new DataFactory[subDataInfo.length];
            for (int i=0; i<factories.length; i++) {
                factories[i] = subDataInfo[i].getDataFactory();
            }
            return factories;
        }
        
        public int getNDataInfo() {
            return subDataInfo.length;
        }
        
        /**
         * Returns the DataInfo corresponding to the group's given wrapped 
         * Data object.
         */
        public DataInfo getSubDataInfo(int i) {
            return subDataInfo[i];
        }
        
        private final DataInfo[] subDataInfo;
    }
}
