package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataInfoFactory;
import etomica.data.DataTag;
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
public class DataGroup implements Data, java.io.Serializable {

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
    
    public static class DataInfoGroup extends DataInfo {
        public DataInfoGroup(String label, Dimension dimension, DataInfo[] subDataInfo) {
            super(label, dimension);
            this.subDataInfo = (DataInfo[])subDataInfo.clone();
        }

        public DataInfoFactory getFactory() {
            return new DataInfoGroupFactory(this);
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
        
        public void addTags(DataTag[] newTags) {
            super.addTags(newTags);
            for (int i=0; i<subDataInfo.length; i++) {
                subDataInfo[i].addTags(newTags);
            }
        }
        
        public void addTag(DataTag newTag) {
            super.addTag(newTag);
            for (int i=0; i<subDataInfo.length; i++) {
                subDataInfo[i].addTag(newTag);
            }
        }
        
        public Data makeData() {
            Data[] subData = new Data[subDataInfo.length];
            for (int i=0; i<subData.length; i++) {
                subData[i] = subDataInfo[i].makeData();
            }
            return new DataGroup(subData);
        }

        protected final DataInfo[] subDataInfo;
    }
    
    public static class DataInfoGroupFactory extends DataInfoFactory {
        protected DataInfoGroupFactory(DataInfoGroup template) {
            super(template);
            subDataInfo = (DataInfo[])template.subDataInfo.clone();
        }
        
        public DataInfo makeDataInfo() {
            DataInfoGroup dataInfo = new DataInfoGroup(label, dimension, subDataInfo);
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
        }
     
        public void setSubDataInfo(DataInfo[] newSubDataInfo) {
            subDataInfo = (DataInfo[])newSubDataInfo.clone();
        }
        
        public DataInfo[] getSubDataInfo() {
            return subDataInfo;
        }
        
        protected DataInfo[] subDataInfo;
    }
}
