package etomica.data;


/**
 * Abstract container of simulation data and information describing it. Data
 * objects may encapsulate any type of data (such as one or more double values),
 * and are transmitted from a DataSource to a DataSink, perhaps passing through
 * a series of data processing elements along the way. The abstract Data class
 * provides methods to obtain descriptive information about the data (its
 * DataInfo), a method to copy the data to another existing data object
 * instance, and a method to make a copy of the data to a new data object
 * instance. Different subclasses of Data are generally not interchangable, and
 * objects processing or receiving Data instances often expect a specific
 * subclass.
 * <p>
 * The data held by a Data instance is mutable, but the structure of a Data
 * instance is not. The "structure" of a Data instance refers to the number and
 * type of data it holds. A DataDoubleArray, for example, holds an array of
 * double; the values in the array may be changed, but the length of the array
 * (the structure) cannot. Thus any object holding a reference to a Data
 * instance can be assured that the number and type of data it holds will not be
 * changed by something else that also references the instance.
 * <p>
 * Information about the data is encapsulated in a DataInfo instance that is
 * held by the Data instance. The DataInfo holds a descriptive label, a
 * Dimension instance that describes the physical dimensions (e.g., mass, time)
 * of the data, and a DataFactory that will make new independent instances of
 * the same Data type and with the same structure. The DataInfo instance is,
 * like the Data structure, completely immutable.
 * 
 * @author David Kofke and Andrew Schultz
 * 
 * @see DataInfo
 */

public interface Data {

    /**
     * Returns a new instance of a data object, formed as a deep copy of this
     * instance.
     */
    public Data makeCopy();

    /**
     * Deep-copies the data from the given object to this one.
     * 
     * @throws ClassCastException
     *             if the given data object is not of a type that can be copied
     *             to this one.
     */
    public void E(Data data);
    
}
