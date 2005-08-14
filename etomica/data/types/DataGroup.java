package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;


/**
 * Gathers one or more Data instances into a single Data object.  The primary
 * use is to enables a data handler to take a single data stream and emit multiple
 * different streams from it.  An example is given by the accumulators, which
 * take data and compile one or more new sets of data from it.
 * </ul>
 * Data sinks typically need to manipulate data groups to extract the individual
 * data elements in them.
 * <p>
 * The DataInfo instance for the DataGroup has a label that describes the collection,
 * and a dimension that is the dimension of all of the data objects it holds (if they
 * are all the same) or UNDEFINED (if they are not).
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
     * Forms a data group from the given array of data objects.  Given data
     * array is cloned, but the data instances it holds are kept by this DataGroup.
     * 
     * @param label a descriptive label for the data collection
     * @param dimension physical dimensions of the data, or UNDEFINED if they are not of the same physical dimensions
     * @param data array of data to be encapsulated in this group
     */
    public DataGroup(String label, Dimension dimension, Data[] data) {
        super(new DataInfo(label, dimension, getFactory(data)));
        this.data = (Data[])data.clone();
    }
    
    /**
     * Copy constructor.  Performs a deep copy, forming a new
     * DataGroup with new instances of copies of this instance's
     * data objects.
     */
    public DataGroup(DataGroup group) {
        super(group);
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
     * Applies the E method to all data elements held, in a one-to-one
     * correspondence with the elements in the given data group. 
     * 
     * @throws ClassCastException if the given data is not an instance of DataGroup.
     * @throws IllegalArgumentException if the given DataGroup has a different number of data elements than this DataGroup.
     * 
     */
    public void E(Data data) {
        if(((DataGroup)data).data.length != this.data.length) {
            throw new IllegalArgumentException("Attempt to copy data groups of different length: (this.length, argument's length): ("+this.data.length+", "+((DataGroup)data).data.length+")");
        }
        for (int i=0; i<this.data.length; i++) {
            this.data[i].E(((DataGroup)data).getData(i));
        }
    }

    /**
     * Returns the i-th data element in the group, counting from 0.
     *
     * @throws ArrayIndexOutOfBounds exception if the given value does not reference a legitimate element
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
        StringBuffer string = new StringBuffer(dataInfo.getLabel());
        for(int i=0; i<data.length; i++) {
            string.append("\n"); //newline?
            string.append(data.toString());
        }
        return string.toString();
    }
    
    private final Data[] data;
    
    /**
     * Returns a factory for a data group holding copies of the
     * given data objects.  Each new instance made by the factory
     * will be a DataGroup having its owns copy of the given data.
     * When the factory is made, the given Data array is cloned, 
     * but the data instances it holds are not.
     */
    public static DataFactory getFactory(Data[] data) {
        return new Factory(data);
    }

    public static class Factory implements DataFactory, Serializable {
        
        final Data[] data;
        
        Factory(Data[] data) {
            this.data = (Data[])data.clone();
        }
        
        public Data makeData(String label, Dimension dimension) {
            Data[] newData = new Data[data.length];
            for(int i=0; i<data.length; i++) {
                if(data[i] != null) {
                    newData[i] = data[i].makeCopy();
                }
            }
            return new DataGroup(label, dimension, newData);
        }
        
        public DataInfo[] getDataInfoArray() {
            DataInfo[] array = new DataInfo[data.length];
            for (int i=0; i<data.length; i++) {
                array[i] = data[i].getDataInfo();
            }
            return array;
        }
        
        public Data[] getData() {
            return (Data[])data.clone();
        }
        
        public Class getDataClass() {
            return DataGroup.class;
        }
        
    }

}
