package etomica.data;

import etomica.DataSource;
import etomica.DataTranslator;
import etomica.exception.MethodNotImplementedException;
import etomica.units.Dimension;
import etomica.utility.Arrays;


/**
 * Collects one or more data sources into a single DataSource. Data returned
 * by the DataSourceGroup is the concatenation of the data returned
 * by the component DataSource instances.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Apr 8, 2005 by kofke
 */
public class DataSourceGroup implements DataSource {

    /**
     * Forms a DataGroup that contains no data sources.  Must
     * populate group with addDataSource method. 
     */
    public DataSourceGroup() {
        super();
        dataSources = new DataSource[0];
        data = new double[0];
        dataLength = new int[0];
    }
    
    /**
     * Forms a DataGroup containing the given data sources.  Given
     * array is copied to another array internally.
     */
    public DataSourceGroup(DataSource[] sources) {
        this();
        if(sources != null) {
            for(int i=0; i<sources.length; i++) {
                addDataSource(sources[i]);
            }
        }
    }

    /**
     * Returns data that is the concatenation of all the data returned
     * by the component groups.  Length of returned array will be equal
     * to the sum of the length of the arrays returned by the component
     * groups at the time this method is called (it is possible for the
     * length to differ from one call to the next, if the length of one
     * or more of the component source's data changes).
     */
    public double[] getData() {
        int cursor = 0;
        for(int i=0; i<dataSources.length; i++) {
            double[] nextData = dataSources[i].getData();
            if(nextData.length != dataLength[i]) {
                adjustDataLength(nextData.length - dataLength[i]);
                dataLength[i] = nextData.length;
                if(!firstCall) System.err.println("Warning: adjusting data length in DataSourceGroup");
            }
            System.arraycopy(nextData, 0, data, cursor, dataLength[i]);
            cursor += dataLength[i];
        }
        firstCall = false;
        return data;
    }
    
    //used by getData
    private void adjustDataLength(int growth) {
        double[] newData = new double[data.length + growth];
        System.arraycopy(data, 0, newData, 0, (growth>0) ? data.length : newData.length);
        data = newData;
    }
    
    /**
     * Adds the given DataSource to those held by the group.  The new data source
     * must give data having the same physical dimensions (as indicated by its
     * getDimension()) as those already added (if any).
     */
    public void addDataSource(DataSource newSource) {
        if(dimension != null && (newSource.getDimension() != dimension)) {
//               || newSource.getTranslator() != translator)) {
            throw new IllegalArgumentException("Illegal attempt to group data sources that give data having different physical dimensions");
        }
        dataSources = (DataSource[])Arrays.addObject(dataSources, newSource);
        dataLength = Arrays.resizeArray(dataLength, dataSources.length);
        dimension = newSource.getDimension();
        translator = newSource.getTranslator();
    }

    /* (non-Javadoc)
     * @see etomica.DataSource#getLabel()
     */
    public String getLabel() {
        if(dataSources.length > 0) return dataSources[0].getLabel();
        else return "";
    }

    /* (non-Javadoc)
     * @see etomica.DataSource#getDimension()
     */
    public Dimension getDimension() {
        return dimension;
    }

    /* (non-Javadoc)
     * @see etomica.DataSource#getTranslator()
     */
    public DataTranslator getTranslator() {
        //FIXME translator needs to separate data pieces and translate each
        throw new MethodNotImplementedException("translator not implemented for DataGroup");
//        return translator;
    }

    private DataSource[] dataSources;
    private Dimension dimension = null;
    private DataTranslator translator = null;
    private double[] data;
    private int[] dataLength;
    private boolean firstCall = true;
}
