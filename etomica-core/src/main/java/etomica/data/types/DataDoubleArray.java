/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import java.util.Arrays;

import etomica.math.function.IFunction;
import etomica.data.DataInfo;
import etomica.data.DataInfoFactory;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataInfoFactory;
import etomica.units.Dimension;

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
//  0     1     2     3     4     5     6     7     8     9    10    11   arrayIndex
//(000) (001) (002) (010) (011) (012) (100) (101) (102) (110) (111) (112) index given to getValue(int[])
//for this example, size = {2, 2, 3}, jumpCount = {6, 3, 1}
//note that number of sites = size[0]*jumpCount[0]

public class DataDoubleArray implements IData, java.io.Serializable {

    /**
     * Constructs a new one-dimensional array of the given length.
     * 
     * @param nValues
     *            length of the one-dimensional array
     */
    public DataDoubleArray(int nValues) {
        this(new int[] {nValues});
    }
    
    /**
     * Constructs a new multidimensional array of the given shape.
     * 
     * @param arrayShape
     *            length of the array in each dimension
     */
    public DataDoubleArray(int[] arrayShape) {
        super();
        this.arrayShape = arrayShape.clone();
        jumpCount = new int[arrayShape.length];
        //row-wise definition, as done in RectangularLattice
        jumpCount[arrayShape.length-1] = 1;
        for(int i=arrayShape.length-1; i>0; i--) {
            jumpCount[i-1] = jumpCount[i]*arrayShape[i];
        }
        x = new double[arrayShape[0]*jumpCount[0]];
    }

   /**
    * Wraps the given array, xData in a DataDoubleArray of the given shape.
    * 
    * @param arrayShape
    *            length of the array in each dimension
    * @param xData
    *            actual data to be held by this instance.  The array is not
    *            cloned.
    */
    public DataDoubleArray(int[] arrayShape, double[] xData) {
        super();
        this.arrayShape = arrayShape.clone();
        jumpCount = new int[arrayShape.length];
        //row-wise definition, as done in RectangularLattice
        jumpCount[arrayShape.length-1] = 1;
        for(int i=arrayShape.length-1; i>0; i--) {
            jumpCount[i-1] = jumpCount[i]*arrayShape[i];
        }
        if (jumpCount[0]*arrayShape[0] != xData.length) {
            throw new IllegalArgumentException("length of xData must be equal to product of arrayShapes");
        }
        x = xData;
    }
    
    /**
     * Returns the dimension of the array, which is the number of integer
     * indices needed to access one of its elements.
     */
    public int getArrayDimension() {
        return jumpCount.length;
    }

    /**
     * Sets the wrapped array of values to the values in the given
     * instance.
     */
    public void E(IData y) {
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
            throw new IllegalArgumentException("DataDoubleArray can only be set equal to another of equal length");
        }
    }

    /**
     * Plus-equals (+=) operation. Element-by-element addition of the values in
     * the given array to those in this one.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void PE(IData y) {
        double[] yx = ((DataDoubleArray) y).x;
        for (int i = 0; i < x.length; i++) {
            x[i] += yx[i];
        }

    }

    /**
     * Minus-equals (-=) operation. Element-by-element subtraction of the values
     * in the given array from those in this one.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void ME(IData y) {
        double[] yx = ((DataDoubleArray) y).x;
        for (int i = 0; i < x.length; i++) {
            x[i] -= yx[i];
        }
    }

    /**
     * Times-equals (*=) operation. Applied element-by-element.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void TE(IData y) {
        double[] yx = ((DataDoubleArray) y).x;
        for (int i = 0; i < x.length; i++) {
            x[i] *= yx[i];
        }

    }

    /**
     * Divide-equals (/=) operation. Applied element-by-element.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if the array in the given object is smaller than this
     *             instance's array.
     */
    public void DE(IData y) {
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
     * Divide-equals (*=) operation. Divides all elements by the given value.
     */
    public void DE(double y) {
        TE(1.0/y);
    }

    /**
     * Replaces all values by the value of the function applied to each.
     */
    public void map(IFunction function) {
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
        return arrayShape[i];
    }

    /**
     * Returns the i-th value of the wrapped array.
     * @throws ArrayIndexOutOfBoundsException if not (0 <= i < getLength()).
     */
    public double getValue(int i) {
        return x[i];
    }

    /**
     * Returns the element of the array indicated by the given set of indices.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if length of given array is less than getArrayDimension (if
     *             length is greater, excess indices are ignored).
     */
    public double getValue(int[] index) {
        return x[arrayIndex(index)];
    }
    
    /**
     * Returns the index in the 1-d array for the value corresponding
     * to the given multidimensional array index.
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
     * Returns the encapsulated array, implementing the Data interface.
     */
    public void assignTo(double[] array) {
        System.arraycopy(x, 0, array, 0, x.length);
    }

    /**
     * Assigns values in the i=th "hypercolumn" of this array to the given array.
     */
    public void assignColumnTo(int i, double[] array) {
        System.arraycopy(x,i*jumpCount[0],array,0,jumpCount.length == 1 ? x.length : jumpCount[0]);
    }

    /**
     * Assigns values in given array to the i-th "hypercolumn" of this data array.
     */
    public void assignColumnFrom(int i, double[] array) {
        System.arraycopy(array,0,x,i*jumpCount[0],jumpCount.length == 1 ? x.length : jumpCount[0]);
    }

    /**
     * Assigns values in the given array to the subsection of this array
     * beginning at the specified index.
     * 
     * @param idx
     *            index array of length less than or equal to getArrayDimension
     *            (right-padded with zeros if shorter)
     * @param array
     *            values to be copied into this array
     */
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
        return x.toString();
    }
    
    private static final long serialVersionUID = 1L;
    private final double[] x;
    private final int[] jumpCount;
    private final int[] arrayShape;
    
    public static class DataInfoDoubleArray extends DataInfo {
        protected final int[] arrayShape;
        
        public DataInfoDoubleArray(String label, Dimension dimension, int[] arrayShape) {
            super(label, dimension);
            this.arrayShape = arrayShape;
        }
        
        /**
         * Returns the size of each dimension in the constructed DataDoubleArray.
         */
        public int[] getArrayShape() {
            return arrayShape.clone();
        }
        
        /**
         * Returns the total number of elements in the constructed DataDoubleArray.
         * This is the product of the elements in the array returned by getArrayShape().
         */
        public int getLength() {
            int n = 1;
            for(int i=arrayShape.length-1; i>=0; i--) {
                n *= arrayShape[i];
            }
            return n;
        }
        
        public IEtomicaDataInfoFactory getFactory() {
            return new DataInfoDoubleArrayFactory(this);
        }

        public IData makeData() {
            return new DataDoubleArray(arrayShape);
        }

        private static final long serialVersionUID = 1L;
    }
    
    public static class DataInfoDoubleArrayFactory extends DataInfoFactory {
        protected DataInfoDoubleArrayFactory(DataInfoDoubleArray template) {
            super(template);
            arrayShape = template.arrayShape.clone();
        }
        
        public IEtomicaDataInfo makeDataInfo() {
            DataInfoDoubleArray dataInfo = new DataInfoDoubleArray(label, dimension, arrayShape);
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
        }
        
        /**
         * Sets the array shape.  The array is copied so further
         * changes made to the given array will not be affect this factory.
         */
        public void setArrayShape(int[] newArrayShape) {
            arrayShape = newArrayShape.clone();
        }
        
        /**
         * Returns a copy of the array shape array.
         */
        public int[] getArrayShape() {
            return arrayShape.clone();
        }
        
        private static final long serialVersionUID = 1L;
        private int[] arrayShape;
    }
}
