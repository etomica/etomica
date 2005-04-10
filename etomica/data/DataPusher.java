package etomica.data;

import etomica.DataSink;
import etomica.DataTranslator;
import etomica.units.Dimension;
import etomica.utility.Arrays;

/**
 * An object that may hold one or more data sinks to which it can 
 * push its data.
 */
public abstract class DataPusher {

    /**
     * Sets the label associated with the data passing through this stream.
     */
    public void setLabel(String label) {
        this.label = label;
    }

    /**
     * Returns the label most recently given to setLabel.  If
     * not previously specified, returns an empty string.  Will
     * not return null.
     */
    public String getLabel() {
        return (label == null) ? "" : label;
    }

    /**
     * Sets label to the given value if it was not previously set.
     * If setLabel was previously called, this method has no effect.
     * This method is usually invoked automatically when this object
     * is placed in a data stream.
     */
    public void setDefaultLabel(String defaultLabel) {
        if(label == null) setLabel(defaultLabel);
    }

    public abstract DataTranslator getTranslator();
    
    /**
     * Method called by subclasses to move data into sinks.
     */
    protected void pushData(double[] data) {
        for(int i=dataSinkList.length-1; i>=0; i--) {
            dataSinkList[i].putData(data);
        }
    }

    /**
     * @return Returns the data sinks.
     */
    public DataSink[] getDataSinks() {
        return dataSinkList;
    }

    /**
     * @param dataSinks The data sinks to set.
     */
    public void setDataSinks(DataSink[] dataSinks) {
        if(dataSinks == null) {
            dataSinkList = new DataSink[0];
            return;
        }
        dataSinkList = (DataSink[])dataSinks.clone();
        for(int i=0; i<dataSinks.length; i++) {
            dataSinkList[i].setDimension(dimension);
            dataSinkList[i].setDefaultLabel(label);
        }
    }

    public void addDataSink(DataSink dataSink) {
        if(dataSink == null) return;
        dataSinkList = (DataSink[])Arrays.addObject(dataSinkList, dataSink);
        dataSink.setDimension(dimension);
        dataSink.setDefaultLabel(label);
    }

    /**
     * Removes the specified data sink from this manager.
     * @param dataSink data sink to be removed from this list, if present.
     * @return <tt>true</tt> if the manager contained the specified data sink.
     */
    public void removeDataSink(DataSink dataSink) {
        dataSinkList = (DataSink[])Arrays.removeObject(dataSinkList, dataSink);
    }
    
    public void setDimension(Dimension dimension) {
        this.dimension = dimension;
        for(int i=0; i<dataSinkList.length; i++) {
            dataSinkList[i].setDimension(dimension);
        }
    }

    public Dimension getDimension() {
        return dimension;
    }

    /**
     * A string describing the property measured by the meter
     */
    protected String label = null;
    private Dimension dimension = Dimension.UNDEFINED;
    protected DataSink[] dataSinkList = new DataSink[0];


}
