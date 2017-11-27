/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import com.fasterxml.jackson.annotation.JsonProperty;
import etomica.data.*;
import etomica.math.function.IFunction;
import etomica.units.dimensions.Dimension;


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
 *
 * @author David Kofke and Andrew Schultz
 */
public class DataGroup implements IData, java.io.Serializable {

    /**
     * Forms a data group from the given array of data objects. Given data array
     * is cloned, but the data instances it holds are kept by this DataGroup.
     * <p>
     * Dimension associated with the DataGroup is Dimension.MIXED if the Data it
     * holds have different dimensions; otherwise dimension is the common
     * dimension of the grouped Data instances.
     * 
     * @param data
     *            array of data to be encapsulated in this group
     * 
     * @throws NullPointerException
     *             if any of the elements of the Data array are null
     */
    public DataGroup(IData[] data) {
        super();
        this.data = data.clone();
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
    public void E(IData newData) {
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
     * @throws ArrayIndexOutOfBoundsException
     *             exception if the given value does not reference a legitimate
     *             element
     */
    public IData getData(int i) {
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
            string.append(data[i].toString());
            if(i < data.length-1) string.append(", ");
        }
        return string.toString();
    }
    
    private static final long serialVersionUID = 1L;
    protected final IData[] data;
    
    public static class DataInfoGroup extends DataInfo {
        public DataInfoGroup(String label, Dimension dimension, IDataInfo[] subDataInfo) {
            super(label, dimension);
            this.subDataInfo = subDataInfo.clone();
        }

        public IDataInfoFactory getFactory() {
            return new DataInfoGroupFactory(this);
        }
        
        public int getLength() {
            int l = 0;
            for (int i=0; i<subDataInfo.length; i++) {
                l += subDataInfo[i].getLength();
            }
            return l;
        }
        
        public int getNDataInfo() {
            return subDataInfo.length;
        }
        
        /**
         * Returns the DataInfo corresponding to the group's given wrapped 
         * Data object.
         */
        public IDataInfo getSubDataInfo(int i) {
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
        
        public IData makeData() {
            IData[] subData = new IData[subDataInfo.length];
            for (int i=0; i<subData.length; i++) {
                subData[i] = subDataInfo[i].makeData();
            }
            return new DataGroup(subData);
        }

        private static final long serialVersionUID = 1L;

        @JsonProperty
        protected final IDataInfo[] subDataInfo;
    }
    
    public static class DataInfoGroupFactory extends DataInfoFactory {

        protected DataInfoGroupFactory(DataInfoGroup template) {
            super(template);
            subDataInfo = template.subDataInfo.clone();
        }
        
        public IDataInfo makeDataInfo() {
            DataInfoGroup dataInfo = new DataInfoGroup(label, dimension, subDataInfo);
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
        }
     
        public void setSubDataInfo(IDataInfo[] newSubDataInfo) {
            subDataInfo = (IDataInfo[])newSubDataInfo.clone();
        }
        
        public IDataInfo[] getSubDataInfo() {
            return subDataInfo;
        }
        
        private static final long serialVersionUID = 1L;
        protected IDataInfo[] subDataInfo;
    }

    public void assignTo(double[] array) {
        int j = 0;
        for (int i=0; i<data.length; i++) {
            for (int k = 0; k < data[i].getLength(); k++) {
                array[j] = data[i].getValue(k);
                j++;
            }
        }
    }

    public void DE(IData y) {
        for (int i=0; i<data.length; i++) {
            data[i].DE(((DataGroup)y).getData(i));
        }
    }

    public void E(double y) {
        for (int i=0; i<data.length; i++) {
            data[i].E(y);
        }
    }

    public int getLength() {
        int l = 0;
        for (int i=0; i<data.length; i++) {
            l += data[i].getLength();
        }
        return l;
    }

    public double getValue(int i) {
        int l = 0;
        for (int j=0; j<data.length; j++) {
            int jl = data[j].getLength();
            l += jl;
            if (l > i) {
                return data[j].getValue(jl-(l-i));
            }
           
        }
        throw new IllegalArgumentException("Length is only "+getLength());
    }

    public boolean isNaN() {
        for (int i=0; i<data.length; i++) {
            if (data[i].isNaN()) {
                return true;
            }
        }
        return false;
    }

    public void map(IFunction function) {
        for (int i=0; i<data.length; i++) {
            data[i].map(function);
        }
    }

    public void ME(IData y) {
        for (int i=0; i<data.length; i++) {
            data[i].ME(((DataGroup)y).getData(i));
        }
    }

    public void PE(IData y) {
        for (int i=0; i<data.length; i++) {
            data[i].PE(((DataGroup)y).getData(i));
        }
    }

    public void PE(double y) {
        for (int i=0; i<data.length; i++) {
            data[i].PE(y);
        }
    }

    public void TE(IData y) {
        for (int i=0; i<data.length; i++) {
            data[i].TE(((DataGroup)y).getData(i));
        }
    }

    public void TE(double y) {
        for (int i=0; i<data.length; i++) {
            data[i].TE(y);
        }
    }

    public void DE(double y) {
        for (int i=0; i<data.length; i++) {
            data[i].DE(y);
        }
    }
}
