package etomica.data;

import java.io.Serializable;

import etomica.Data;
import etomica.DataInfo;
import etomica.units.Dimension;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 16, 2005 by kofke
 */
public class DataGroup extends Data {

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


    public void E(Data data) {
        for (int i=0; i<this.data.length; i++) {
            this.data[i].E(((DataGroup)data).getData(i));
        }
    }

    public Data getData(int i) {
        return data[i];
    }
    
    public int getNData() {
        return data.length;
    }
    
    public String toString() {
        StringBuffer string = new StringBuffer(dataInfo.getLabel());
        for(int i=0; i<data.length; i++) {
            string.append("\n"); //newline?
            string.append(data.toString());
        }
        return string.toString();
    }
    
    final Data[] data;
    
    public static DataFactory getFactory(Data[] data) {
        return new Factory(data);
    }

    private static class Factory implements DataFactory, Serializable {
        
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
        
        public Class getDataClass() {
            return DataGroup.class;
        }
        
        public DataFactory copy() {
            return DataGroup.getFactory(data);
        }
        
    }

}
