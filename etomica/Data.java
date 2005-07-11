package etomica;


/**
 * Abstract container of simulation data and information describing it. Data
 * objects may encapsulate any type of data (such as one or more double values),
 * and proceed from a DataSource to a DataSink, perhaps passing through a series
 * of data processing elements along the way. The abstract Data class provides
 * methods to obtain descriptive information about the data (its physical
 * dimensions and a descriptive label), a method to copy the data to another
 * existing data object instance, and a method to make a copy of the data to a
 * new data object instance. Different subclasses of Data are generally not
 * interchangable, and objects processing or receiving Data instances often
 * expect a specific subclass.
 * 
 * @author David Kofke and Andrew Schultz
 * 
 * @see DataInfo
 */

/*
 * History Created on Jun 15, 2005 by kofke
 */
public abstract class Data implements java.io.Serializable {

    /**
     * Constructs the object with the given DataInfo instance, which is returned
     * by the getDataInfo method. DataInfo is declared final and cannot be
     * changed after construction.
     */
    public Data(DataInfo dataInfo) {
        this.dataInfo = dataInfo;
    }

    /**
     * Copy constructor, used by subclasses.  Constructs instance by making
     * a copy of the DataInfo in the given Data instance.
     */
    protected Data(Data data) {
        this(new DataInfo(data.dataInfo));
    }

    /**
     * Returns a DataInfo object, which contains descriptive information about
     * the data held by this object.
     */
    public DataInfo getDataInfo() {
        return dataInfo;
    }

    /**
     * Returns a new instance of a data object, formed as a deep copy of this
     * instance.
     */
    public abstract Data makeCopy();

    /**
     * Deep-copies the data from the given object to this one.
     * 
     * @throws ClassCastException
     *             if the given data object is not of a type that can be copied
     *             to this one.
     */
    public abstract void E(Data data);

    protected final DataInfo dataInfo;
    
 }
