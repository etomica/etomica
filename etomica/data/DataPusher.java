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

    public void setLabel(String label) {
        this.label = label;
    }

    public String getLabel() {
        return label;
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
//            if(dataSinkList[i].getDimension() == Dimension.NULL) {
                dataSinkList[i].setDimension(dimension);
//            }
//            if(dataSinkList[i].getLabel().equals("")) {
                dataSinkList[i].setLabel(label);
//            }
        }
    }

    public void addDataSink(DataSink dataSink) {
        if(dataSink == null) return;
        dataSinkList = (DataSink[])Arrays.addObject(dataSinkList, dataSink);
        dataSink.setDimension(dimension);
        dataSink.setLabel(label);
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
    protected String label = "";
    private Dimension dimension = Dimension.UNDEFINED;
    protected DataSink[] dataSinkList = new DataSink[0];


}
