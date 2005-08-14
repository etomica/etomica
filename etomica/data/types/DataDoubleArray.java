package etomica.data.types;

import java.io.Serializable;
import java.util.Arrays;

import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.units.Dimension;
import etomica.utility.Function;

/**
 * Data object that wraps a mutable array of doubles. Data may be treated alternately as a simple 
 * one-dimensional array accessed via the getData() method, or as a multidimensional array
 * with indexed elements accessed via the getData(int[]) method.  Shape of multidimensional array
 * is specified at construction and cannot be changed afterwards. 
 * <p>
 * All arithmetic methods throw ClassCastException if given a Data instance that
 * is not a DataDoubleArray, and throw an IllegalArgumentException if the given DataDoubleArray is
 * not of the same length as this one (array shape is not considered for arithmetic operations; only
 * total length is relevant).
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

/*
 * History Created on Jun 15, 2005
 */

//Internally, sites are stored in a 1-D array of objects, and are accessed by
//unrolling the index specification to determine the storage-array index.

//Example showing internal ordering of elements
//0     1     2     3     4     5     6     7     8     9    10    11   arrayIndex
//(000) (001) (002) (010) (011) (012) (100) (101) (102) (110) (111) (112) index given to getValue(int[])
//for this example, size = {3, 2, 2}, jumpCount = {6, 3, 1}
//note that number of sites = size[n]*jumpCount[n] where n = size.length-1

public class DataDoubleArray extends Data implements DataArithmetic {

    /**
     * Constructs a new instance with the given DataInfo.
     */
    public DataDoubleArray(String label, Dimension dimension, int nValues) {
        this(label, dimension, new int[] {nValues});
    }
    
    
    public DataDoubleArray(String label, Dimension dimension, int[] arrayShape) {
        super(new DataInfo(label, dimension, getFactory(arrayShape)));
        jumpCount = (int[])arrayShape.clone();
        //row-wise definition, as done in RectangularLattice
        jumpCount[arrayShape.length-1] = 1;
        for(int i=arrayShape.length-1; i>0; i--) {
            jumpCount[i-1] = jumpCount[i]*arrayShape[i];
        }
        x = new double[arrayShape[0]*jumpCount[0]];
    }

    /**
     * Constructs a new instance using the label and dimension from the given DataInfo.
     * New instance will have a new DataInfo; the factory in the new DataInfo will construct
     * DataDoubleArray instances using the given arrayShape (which likely differs from the
     * arrayShape made by the factory in the given DataInfo).
     * <p>
     * This is a convenience constructor, and is equivalent to<br><code>
     * this(dataInfo.getLabel(), dataInfo.getDimension(), arrayShape);</code>
     *  
     * @param dataInfo provides label and Dimension for new instance
     * @param arrayShape specifies dimensions of new data array.
     */
    public DataDoubleArray(DataInfo dataInfo, int[] arrayShape) {
        this(dataInfo.getLabel(), dataInfo.getDimension(), arrayShape);
    }

    /**
     * Copy constructor.  Makes a new DataDoubleArray having the same DataInfo instance
     * as the given DataDoubleArray.  New instance will encapsulate a double array that
     * is a clone of the one in the given instance. 
     */
    public DataDoubleArray(DataDoubleArray data) {
        super(data);
        x = (double[])data.x.clone();
        jumpCount = data.jumpCount;
    }

    /**
     * Returns a copy of this instance. Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataDoubleArray(this);
    }
    
    public int getArrayDimension() {
        return jumpCount.length;
    }

    /**
     * Sets the wrapped array of values to the values in the given
     * instance.
     */
    public void E(Data y) {
        this.E(((DataDoubleArray) y).x);
    }

    /**
     * Copies each value in the given array to the values in the wrapped
     * array.
     * 
     * @throws IllegalArgumentException if the arrays are of different lengths
     */
    public void E(double[] y) {
        if (y.length == x.length) {
            System.arraycopy(y, 0, x, 0, x.length);
        } else {
            throw new IllegalArgumentException("DataDoubleArray can only be set equal to another of equal arrayShape");
        }
    }

    /**
     * Plus-equals (+=) operation. Element-by-element addition of the values in
     * the given array to those in this one.
     * 
     * @throws ArrayIndexOutOfBoundException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void PE(DataArithmetic y) {
        double[] yx = ((DataDoubleArray) y).x;
        for (int i = 0; i < x.length; i++) {
            x[i] += yx[i];
        }

    }

    /**
     * Minus-equals (-=) operation. Element-by-element subtraction of the values
     * in the given array from those in this one.
     * 
     * @throws ArrayIndexOutOfBoundException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void ME(DataArithmetic y) {
        double[] yx = ((DataDoubleArray) y).x;
        for (int i = 0; i < x.length; i++) {
            x[i] -= yx[i];
        }
    }

    /**
     * Times-equals (*=) operation. Applied element-by-element.
     * 
     * @throws ArrayIndexOutOfBoundException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void TE(DataArithmetic y) {
        double[] yx = ((DataDoubleArray) y).x;
        for (int i = 0; i < x.length; i++) {
            x[i] *= yx[i];
        }

    }

    /**
     * Divide-equals (/=) operation. Applied element-by-element.
     * 
     * @throws ArrayIndexOutOfBoundException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void DE(DataArithmetic y) {
        double[] yx = ((DataDoubleArray) y).x;
        for (int i = 0; i < x.length; i++) {
            x[i] /= yx[i];
        }

    }

    /**
     * Equals (=) operation. Sets all elements to given value.
     */
    public void E(double y) {
        Arrays.fill(x, y);
    }

    /**
     * Plus-equals (+=) operation. Increments all elements by the given value.
     */
    public void PE(double y) {
        for (int i = 0; i < x.length; i++) {
            x[i] += y;
        }
    }

    /**
     * Times-equals (*=) operation. Multiplies all elements by the given value.
     */
    public void TE(double y) {
        for (int i = 0; i < x.length; i++) {
            x[i] *= y;
        }
    }

    /**
     * Replaces all values by the value of the function applied to each.
     */
    public void map(Function function) {
        for (int i = 0; i < x.length; i++) {
            x[i] = function.f(x[i]);
        }
    }

    /**
     * Returns the total length of the wrapped array.  This is the
     * product of the lengths in each dimension, and is the length
     * of the wrapped 1D array that stores the array values.
     */
    public int getLength() {
        return x.length;
    }
    
    /**
     * Returns the length of the array in the i-th dimension.
     */
    public int getArrayShape(int i) {
        return ((Factory)dataInfo.getDataFactory()).arrayShape[i];
    }

    /**
     * Returns the i-th value of the wrapped array.
     * @throws ArrayIndexOutOfBoundsException if not (0 <= i < getLength()).
     */
    public double getValue(int i) {
        return x[i];
    }

    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#site(int[])
     */
    public double getValue(int[] index) {
        return x[arrayIndex(index)];
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
     * Returns the encapsulated array.
     */
    public double[] getData() {
        return x;
    }

    /**
     * Returns true if any element of the encapsulated array is not-a-number, as
     * given by Double.isNaN.
     */
    public boolean isNaN() {
        for (int i = 0; i < x.length; i++) {
            if (Double.isNaN(x[i]))
                return true;
        }
        return false;
    }

    /**
     * Returns the encapsulated array, implementing the DataArithmetic interface.
     */
    public void assignTo(double[] array) {
        System.arraycopy(x, 0, array, 0, x.length);
    }

    public void assignColumnTo(int i, double[] array) {
        System.arraycopy(x,i*jumpCount[0],array,0,jumpCount.length == 1 ? x.length : jumpCount[0]);
    }
    
    public void assignColumnFrom(int i, double[] array) {
        System.arraycopy(array,0,x,i*jumpCount[0],jumpCount.length == 1 ? x.length : jumpCount[0]);
    }

    public void assignSubsectionFrom(int[] idx, double[] array) {
        int offset = 0;
        for (int i=0; i<idx.length; i++) {
            offset += idx[i]*jumpCount[i];
        }
        System.arraycopy(array,0,x,offset,array.length);
    }

    /**
     * Returns a string formed from the dataInfo label and the array values.
     */
    public String toString() {
        return dataInfo.getLabel() + " " + x.toString();
    }
    
    private final double[] x;
    private final int[] jumpCount;
    //note that arrayShape is available via dataInfo.getDataFactory
    
    public static DataFactory getFactory(int[] arrayShape) {
        return new Factory(arrayShape);
    }

    public static class Factory implements DataFactory, Serializable {
        
        private final int[] arrayShape;
        
        Factory(int[] arrayShape) {
            this.arrayShape = arrayShape;
        }
        
        public Data makeData(String label, Dimension dimension) {
            DataDoubleArray data = new DataDoubleArray(label, dimension, arrayShape);
            return data;
        }
        
        public Class getDataClass() {
            return DataDoubleArray.class;
        }
        
        /**
         * Returns the size of each dimension in the constructed DataDoubleArray.
         */
        public int[] getArrayShape() {
            return (int[])arrayShape.clone();
        }
        
        /**
         * Returns the total number of elements in the constructed DataDoubleArray.
         * This is the product of the elements in the array returned by getArrayShape().
         */
        public int getArrayLength() {
            int n = 1;
            for(int i=arrayShape.length-1; i>=0; i--) {
                n *= arrayShape[i];
            }
            return n;
        }
        
    }

}