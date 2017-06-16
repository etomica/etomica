/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.HashSet;

import etomica.math.function.IFunction;
import etomica.data.DataInfo;
import etomica.data.DataInfoFactory;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataInfoFactory;
import etomica.units.Dimension;

/**
 * Data object that wraps a array of BigDecimals.  Precision is specified at
 * construction and operations will maintain that level of precision.  Data may
 * be treated alternately as a simple one-dimensional array accessed via the
 * getData() method, or as a multidimensional array with indexed elements
 * accessed via the getData(int[]) method.  Shape of multidimensional array
 * is specified at construction and cannot be changed afterwards. 
 * <p>
 * All arithmetic methods throw ClassCastException if given a Data instance
 * that is not a DataDoubleDBArray, and throw an IllegalArgumentException if
 * the given DataDoubleDBArray is not of the same length as this one (array
 * shape is not considered for arithmetic operations; only total length is
 * relevant).
 * 
 * @author David Kofke and Andrew Schultz
 */
public class DataDoubleBDArray implements IData, java.io.Serializable {

    /**
     * Constructs a new one-dimensional array of the given length.
     * 
     * @param nValues
     *            length of the one-dimensional array
     */
    public DataDoubleBDArray(int nValues, int precision) {
        this(new int[] {nValues}, precision);
    }
    
    /**
     * Constructs a new multidimensional array of the given shape.
     * 
     * @param arrayShape
     *            length of the array in each dimension
     */
    public DataDoubleBDArray(int[] arrayShape, int precision) {
        super();
        this.arrayShape = arrayShape.clone();
        jumpCount = new int[arrayShape.length];
        //row-wise definition, as done in RectangularLattice
        jumpCount[arrayShape.length-1] = 1;
        for(int i=arrayShape.length-1; i>0; i--) {
            jumpCount[i-1] = jumpCount[i]*arrayShape[i];
        }
        mc = new MathContext(precision);
        x = new BigDecimal[arrayShape[0]*jumpCount[0]];
        BigDecimal zero = new BigDecimal(0);
        for (int i=0; i<x.length; i++) {
            x[i] = zero;
        }
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
        this.E(((DataDoubleBDArray) y).x);
    }

    /**
     * Copies each value in the given array to the values in the wrapped
     * array.
     * 
     * @throws IllegalArgumentException if the arrays are of different lengths
     */
    public void E(BigDecimal[] y) {
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
        BigDecimal[] yx = ((DataDoubleBDArray) y).x;
        for (int i = 0; i < x.length; i++) {
            if (zeros.contains(yx[i])) continue;
            x[i] = x[i].add(yx[i], mc);
        }
    }
    
    public static HashSet<BigDecimal> zeros = new HashSet<BigDecimal>();
    public static HashSet<BigDecimal> ones = new HashSet<BigDecimal>();
    
    public static void addBDZero(BigDecimal zero) {
        if (!zeros.contains(zero)) zeros.add(zero);
    }

    public static void addBDOne(BigDecimal one) {
        if (!ones.contains(one)) ones.add(one);
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
        BigDecimal[] yx = ((DataDoubleBDArray) y).x;
        for (int i = 0; i < x.length; i++) {
            if (zeros.contains(yx[i])) continue;
            x[i] = x[i].subtract(yx[i], mc);
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
        BigDecimal[] yx = ((DataDoubleBDArray) y).x;
        for (int i = 0; i < x.length; i++) {
            if (ones.contains(yx[i]) || (zeros.contains(yx[i]) && zeros.contains(x[i]))) continue;
            x[i] = x[i].multiply(yx[i], mc);
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
        BigDecimal[] yx = ((DataDoubleBDArray) y).x;
        for (int i = 0; i < x.length; i++) {
            // we can not represent NaN or infinity
            if (ones.contains(yx[i]) || zeros.contains(yx[i]) || yx[i].doubleValue() == 0) continue;
            x[i] = x[i].divide(yx[i], mc);
        }
    }

    /**
     * Equals (=) operation. Sets all elements to given value.
     */
    public void E(double y) {
        if (Double.isNaN(y)) {
            // pretend like it didn't happen
            // we probably initializing stuff to NaN
            return;
        }
        if (Double.isInfinite(y)) {
            throw new RuntimeException("infinity not allowed");
        }
        BigDecimal ybd = new BigDecimal(y, mc);
        Arrays.fill(x, ybd);
    }

    /**
     * Plus-equals (+=) operation. Increments all elements by the given value.
     */
    public void PE(double y) {
        if (y==0) return;
        BigDecimal ybd = new BigDecimal(y, mc);
        for (int i = 0; i < x.length; i++) {
            x[i] = x[i].add(ybd, mc);
        }
    }

    /**
     * Times-equals (*=) operation. Multiplies all elements by the given value.
     */
    public void TE(double y) {
        if (y==1) return;
        BigDecimal ybd = new BigDecimal(y, mc);
        for (int i = 0; i < x.length; i++) {
            x[i] = x[i].multiply(ybd, mc);
        }
    }

    /**
     * Times-equals (*=) operation. Multiplies all elements by the given value.
     */
    public void DE(double y) {
        if (y==1) return;
        BigDecimal ybd = new BigDecimal(y, mc);
        for (int i = 0; i < x.length; i++) {
            x[i] = x[i].divide(ybd, mc);
        }
    }

    /**
     * Replaces all values by the value of the function applied to each.
     */
    public void map(IFunction function) {
        for (int i = 0; i < x.length; i++) {
            x[i] = new BigDecimal(function.f(x[i].doubleValue()), mc);
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
        return x[i].doubleValue();
    }

    /**
     * Returns the element of the array indicated by the given set of indices.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if length of given array is less than getArrayDimension (if
     *             length is greater, excess indices are ignored).
     */
    public double getValue(int[] index) {
        return x[arrayIndex(index)].doubleValue();
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
    public BigDecimal[] getData() {
        return x;
    }

    public boolean isNaN() {
        // BigDecimal may not represent NaN
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
    public void assignColumnTo(int i, BigDecimal[] array) {
        System.arraycopy(x,i*jumpCount[0],array,0,jumpCount.length == 1 ? x.length : jumpCount[0]);
    }

    /**
     * Assigns values in given array to the i-th "hypercolumn" of this data array.
     */
    public void assignColumnFrom(int i, BigDecimal[] array) {
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
    private final BigDecimal[] x;
    private final int[] jumpCount;
    private final int[] arrayShape;
    private final MathContext mc;
    
    public static class DataInfoDoubleBDArray extends DataInfo {
        protected final int[] arrayShape;
        protected final int precision;
        
        public DataInfoDoubleBDArray(String label, Dimension dimension, int[] arrayShape, int precision) {
            super(label, dimension);
            this.arrayShape = arrayShape;
            this.precision = precision;
        }
        
        /**
         * Returns the size of each dimension in the constructed DataDoubleArray.
         */
        public int[] getArrayShape() {
            return arrayShape.clone();
        }
        
        public int getPrecision() {
            return precision;
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
            return new DataInfoDoubleBDArrayFactory(this);
        }

        public IData makeData() {
            return new DataDoubleBDArray(arrayShape, precision);
        }

        private static final long serialVersionUID = 1L;
    }
    
    public static class DataInfoDoubleBDArrayFactory extends DataInfoFactory {
        protected DataInfoDoubleBDArrayFactory(DataInfoDoubleBDArray template) {
            super(template);
            arrayShape = template.arrayShape.clone();
            precision = template.precision;
        }
        
        public IEtomicaDataInfo makeDataInfo() {
            DataInfoDoubleBDArray dataInfo = new DataInfoDoubleBDArray(label, dimension, arrayShape, precision);
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
        
        public void setPrecision(int newPrecision) {
            precision = newPrecision;
        }
        
        public int getPrecision() {
            return precision;
        }
        
        private static final long serialVersionUID = 1L;
        private int[] arrayShape;
        private int precision;
    }
}
