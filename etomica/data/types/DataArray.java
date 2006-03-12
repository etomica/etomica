package etomica.data.types;

import java.io.Serializable;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;

/**
 * Data object that wraps an array of similarly-formed Data instances. All
 * Data instances are formed from the same DataFactory, and differ only in the
 * data they hold.  All hold the same instance of DataInfo.
 * <p>
 * Array may be treated alternately as a simple one-dimensional array accessed
 * via the getData() method, or as a multidimensional array with indexed
 * elements accessed via the getData(int[]) method. Shape of multidimensional
 * array is specified at construction and cannot be changed afterwards.
 * <p>
 * This class might be used, for example, to represent a vector field as a function
 * of x,y,z coordinates.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 25, 2005
 */

//Array structure is the same as that used by DataDoubleArray.
//Internally, sites are stored in a 1-D array of objects, and are accessed by
//unrolling the index specification to determine the storage-array index.

//Example showing internal ordering of elements
//0     1     2     3     4     5     6     7     8     9    10    11   arrayIndex
//(000) (001) (002) (010) (011) (012) (100) (101) (102) (110) (111) (112) index given to getValue(int[])
//for this example, size = {2, 2, 3}, jumpCount = {6, 3, 1}
//note that number of sites = size[0]*jumpCount[0]

public class DataArray extends Data {

    /**
     * Constructs a 1-dimensional array of Data.
     * 
     * @param label
     *            descriptive label used by the common DataInfo of all arrayed
     *            Data instances
     * @param dimension
     *            physical dimensions of the Data held by each array element
     * @param nValues
     *            length of the 1-D array of Data
     * @param arrayElementFactory
     *            constructs a prototype Data element that is copied to populate
     *            the array
     */
    public DataArray(String label, Dimension dimension, int nValues, DataFactory arrayElementFactory) {
        this(label, dimension, new int[] {nValues}, arrayElementFactory);
    }
    
    /**
     * Constructs an arbitrarily shaped array of Data.
     * 
     * @param label
     *            descriptive label used by the common DataInfo of all arrayed
     *            Data instances
     * @param dimension
     *            physical dimensions of the Data held by each array element
     * @param arrayShape
     *            dimensions of the Data array
     * @param arrayElementFactory
     *            constructs a prototype Data element that is copied to populate
     *            the array
     */
    public DataArray(String label, Dimension dimension, int[] arrayShape, DataFactory arrayElementFactory) {
        this(label, dimension, (Factory)getFactory(arrayShape, arrayElementFactory));
    }
    
    /*
     * Used by factory and constructor above
     */
    private DataArray(String label, Dimension dimension, Factory factory) {
        super(new DataInfo(label+" Array", dimension, factory));
        jumpCount = (int[])factory.arrayShape.clone();
        //row-wise definition, as done in RectangularLattice
        if(factory.arrayShape.length > 0) {
            jumpCount[factory.arrayShape.length-1] = 1;
            for(int i=factory.arrayShape.length-1; i>0; i--) {
                jumpCount[i-1] = jumpCount[i]*factory.arrayShape[i];
            }
            dataArray = new Data[factory.arrayShape[0]*jumpCount[0]];
        } else {
            dataArray = new Data[0];
        }
        
        //all array elements get the same DataInfo instance
        Data prototype = factory.arrayElementFactory.makeData(label, dimension);
        for(int i=0; i<dataArray.length; i++) {
            dataArray[i] = prototype.makeCopy();
        }
    }

//    /**
//     * Constructs a new instance using the label and dimension from the given DataInfo.
//     * New instance will have a new DataInfo; the factory in the new DataInfo will construct
//     * DataArray instances using the given arraySize (which likely differs from the
//     * arraySize made by the factory in the given DataInfo).
//     * <p>
//     * This is a convenience constructor, and is equivalent to<br><code>
//     * this(dataInfo.getLabel(), dataInfo.getDimension(), arraySize, ((Factory)dataInfo.getDataFactory()).getArrayElementFactory()));</code>
//     *  
//     * @param dataInfo provides label, Dimension, and arrayElementFactory for new instance
//     * @param arraySize specifies dimensions of new data array.
//     */
    //unclear which factory to use here; assume DataInfo is for a DataArray and use its elementFactory?
    //or take dataInfo as for an arbitrary Data instance, and use its factory?
//    public DataArray(DataInfo dataInfo, int[] arraySize) {
//        this(dataInfo.getLabel(), dataInfo.getDimension(), arraySize,
//                ((Factory)dataInfo.getDataFactory()).getArrayElementFactory());
//    }

    /**
     * Copy constructor.  Makes a new DataArray having the same DataInfo instance
     * as the given DataArray.  New instance will encapsulate a Data array that
     * is a clone of the one in the given instance, with all Data elements made as
     * copies of the original ones (and thus having the same DataInfo). 
     */
    public DataArray(DataArray data) {
        super(data);
        dataArray = new Data[data.dataArray.length];
        for(int i=0; i<dataArray.length; i++) {
            dataArray[i] = data.dataArray[i].makeCopy();
        }
        jumpCount = data.jumpCount;
    }

    /**
     * Returns a copy of this instance. Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataArray(this);
    }
    
    /**
     * Returns the Class of the Data elements held in the array.
     */
    public Class getArrayElementType() {
        return ((Factory)dataInfo.getDataFactory()).getArrayElementFactory().getClass();
    }
    
    /**
     * Returns the dimension of the array, which is the number of integer indices needed
     * to access one of its elements. 
     */
    public int getArrayDimension() {
        return jumpCount.length;
    }

    /**
     * Sets the wrapped array of values to the values in the given
     * instance.
     */
    public void E(Data y) {
        for(int i=0; i<dataArray.length; i++) {
            this.dataArray[i].E(((DataArray)y).dataArray[i]);
        }
    }

    
    /**
     * Returns the length of the array in the i-th dimension.
     */
    public int getArraySize(int i) {
        return ((Factory)dataInfo.getDataFactory()).arrayShape[i];
    }

    /**
     * Returns the i-th value of the wrapped array.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if not (0 <= i < getLength()).
     */
    public Data getData(int i) {
        return dataArray[i];
    }
    
    /**
     * Returns the length of the wrapped array of Data.
     */
    public int getLength() {
        return dataArray.length;
    }

    /**
     * Returns the element of the array indicated by the given set of indices.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if length of given array is less than getArrayDimension (if
     *             length is greater, excess indices are ignored).
     */
    public Data getData(int[] index) {
        return dataArray[arrayIndex(index)];
    }
    
    /**
     * Returns the index in the 1-d array for the site corresponding
     * to the given lattice index.
     */
    private final int arrayIndex(int[] index) {
        int idx = 0;
        for(int i=0; i<jumpCount.length; i++) {
            idx += index[i]*jumpCount[i];
        }
        return idx;
    }

    /**
     * Returns a string formed from the dataInfo label and the array values.
     */
    public String toString() {
        return dataInfo.getLabel() + " " + dataArray.toString();
    }
    
    private final Data[] dataArray;
    private final int[] jumpCount;
    
    /**
     * Returns a new DataFactory that makes DataArray instances of the given shape with
     * elements formed using the given DataFactory.
     */
    public static DataFactory getFactory(int[] arrayShape, DataFactory factory) {
        return new Factory(arrayShape, factory);
    }

    /**
     * DataFactory the constructs DataArray instances of a specific shape.
     * Instantiate using the static DataArray.getFactory method.
     */
    public static class Factory implements DataFactory, Serializable {
        
        protected final int[] arrayShape;
        protected final DataFactory arrayElementFactory;
        
        Factory(int[] arrayShape, DataFactory arrayElementFactory) {
            this.arrayShape= (int[])arrayShape.clone();
            this.arrayElementFactory = arrayElementFactory;
        }
        
        /**
         * Constructs a new DataArray, using the given label and dimension for its
         * DataInfo, and with shape and element type as set for this factory at
         * construction.
         */
        public Data makeData(String label, Dimension dimension) {
            DataArray data = new DataArray(label, dimension, arrayShape, arrayElementFactory);
            return data;
        }
        
        /**
         * Returns DataArray.class, indicating that this factory makes DataArray instances.
         */
        public Class getDataClass() {
            return DataArray.class;
        }
        
        /**
         * Returns the size of each dimension in the constructed DataDoubleArray.
         */
        public int[] getArrayShape() {
            return (int[])arrayShape.clone();
        }
        
        /**
         * Returns the total number of elements in the constructed DataDoubleArray.
         * This is the product of the elements in the array returned by getArraySize().
         */
        public int getArrayLength() {
            if(arrayShape.length == 0) {
                return 0;
            }
            int n = 1;
            for(int i=arrayShape.length-1; i>=0; i--) {
                n *= arrayShape[i];
            }
            return n;
        }
        
        /**
         * Returns the factory that makes the Data object held in the array.
         */
        public DataFactory getArrayElementFactory() {
            return arrayElementFactory;
        }
        
    }//end of Factory

}