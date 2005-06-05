package etomica.data;

import etomica.DataSink;
import etomica.DataSource;
import etomica.DataTranslator;
import etomica.units.Dimension;

/**
 * DataSink that receives and holds data without pushing it downstream. May
 * notify a DataBinManager when the data is changed, which can then access the
 * data via the DataSource methods implemented by this class. Dimension (e.g.,
 * length, time) is set at construction and cannot be changed (this is because
 * the class does not push its data, but instead serves as a data source to its
 * DataBinManager; it therefore does not have a way to notify downstream sinks
 * if its dimension field changes).
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Apr 9, 2005 by kofke
 */
public class DataBin implements DataSink, DataSource {

    /**
     * Constructs DataBin with null DataBinManager
     * 
     * @param dimension
     *            The physical dimensions of the data in this bin.
     */
    public DataBin(Dimension dimension) {
        this(null, dimension);
    }

    /**
     * Construct new DataBin that will notify the given DataBinManager (via the
     * manager's dataChangeNotify method) (if not null) any time the bin's
     * putData method is invoked.
     * 
     * @param dataBinManager
     *            manger of this bin; declared final; null value is permitted
     * @param dimension
     *            The physical dimensions (e.g., length, time) of the data in
     *            this bin. Used to select units (e.g., Anstroms, picoseconds)
     *            when data is displayed or written to file. Declared final.
     */
    public DataBin(DataBinManager dataBinManager, Dimension dimension) {
        y = new double[0];
        this.dimension = dimension;
        this.dataBinManager = dataBinManager;
    }

    /**
     * Stores the given array of values. Data are held in a copy array, so any
     * subsequent changes to given array will not affect values held by the bin.
     */
    public void putData(double[] values) {
        if (y.length != values.length) {
            y = (double[]) values.clone();
        } else {
            System.arraycopy(values, 0, y, 0, values.length);
        }
        if (dataBinManager != null) {
            dataChanged = true;
            dataBinManager.dataChangeNotify(this);
        }
    }

    /**
     * Returns the array holding the data most recently given to putData.
     * Returns a zero-length array if putData was not previously invoked.
     */
    public double[] getData() {
        return y;
    }

    /**
     * Returns the length of the data, as last given to putData.
     */
    public int getDataLength() {
        return y.length;
    }

    /**
     * Sets the physical dimensions (e.g., length, time) of the data in this
     * bin. This is set at construction, so calls to this method are not needed.
     * Exception is thrown if given dimension does not match dimension set in
     * constructor.
     */
    public void setDimension(Dimension dimension) {
        if (this.dimension != dimension)
            throw new IllegalArgumentException(
                    "Cannot change dimension of DataBin after construction");
    }

    /**
     * Returns the physical dimensions (e.g., length, time) of the data in this
     * bin. Used to select units (e.g., Anstroms, picoseconds) when data is
     * displayed or written to file.
     */
    public Dimension getDimension() {
        return dimension;
    }

    /**
     * Sets label to the given value if it was not previously set. If setLabel
     * was previously called, this method has no effect.
     */
    public void setDefaultLabel(String defaultLabel) {
        if (label == null)
            setLabel(defaultLabel);
    }

    /**
     * Defines a label to associate with the data in this bin. Used when data is
     * displayed or written to file.
     */
    public void setLabel(String label) {
        this.label = label;
    }

    /**
     * Returns the label associated with this data. Used when data is displayed
     * or written to file.
     */
    public String getLabel() {
        return (label == null) ? "" : label;
    }

    /**
     * @return Returns the dataTranslator.
     */
    public DataTranslator getTranslator() {
        return dataTranslator;
    }

    /**
     * @param dataTranslator
     *            The dataTranslator to set.
     */
    public void setTranslator(DataTranslator dataTranslator) {
        this.dataTranslator = dataTranslator;
    }

    private double[] y;
    private final Dimension dimension;
    private String label = null;
    private DataTranslator dataTranslator;
    private final DataBinManager dataBinManager;

    //flag used by DataBinManager to keep track of bin changes
    boolean dataChanged = true;
}