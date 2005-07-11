package etomica.data.types;

import java.io.Serializable;
import java.util.Arrays;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.DataFactory;
import etomica.units.Dimension;
import etomica.utility.Function;

/**
 * Data object that wraps a mutable array of doubles, which may be accessed via
 * the getData method.<br>
 * All arithmetic methods throw ClassCastException if given a Data instance that
 * is not of this type.
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

/*
 * History Created on Jun 15, 2005
 */
public class DataDoubleArray extends Data implements DataArithmetic {

    /**
     * Constructs a new instance with the given DataInfo.
     */
    public DataDoubleArray(String label, Dimension dimension) {
        super(new DataInfo(label, dimension, getFactory(0)));
    }

    /**
     * Copy constructor.
     */
    public DataDoubleArray(DataDoubleArray data) {
        super(data);
        int n = data.x.length;
        ((Factory)dataInfo.getDataFactory()).dataLength = n;
        x = new double[n];
        System.arraycopy(data.x, 0, x, 0, x.length);
    }

    /**
     * Returns a copy of this instance. Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataDoubleArray(this);
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
     * array. Adjusts size of this array to that of the given array, if needed.
     */
    public void E(double[] y) {
        if (y.length == x.length) {
            System.arraycopy(y, 0, x, 0, x.length);
        } else {
            x = (double[]) y.clone();
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
     * Replaces wrapped array with a new one of the given length,
     * discarding current values, if parameter is not equal to current length of
     * array.
     * 
     * @throws IllegalArgumentException
     *             if n < 0
     */
    public void setLength(int n) {
        if (n < 0)
            throw new IllegalArgumentException("Illegal length: " + n);
        if (n != x.length)
            ((Factory)dataInfo.getDataFactory()).dataLength = n;
            x = new double[n];
    }

    /**
     * Returns the length of the wrapped array.
     */
    public int getLength() {
        return x.length;
    }

    /**
     * Returns the i-th value of the wrapped array.
     * @throws ArrayIndexOutOfBoundsException if not (0 <= i < getLength()).
     */
    public double getValue(int i) {
        return x[i];
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
    public double[] toArray() {
        return x;
    }

    public DataArithmetic toArithmetic(DataArithmetic data) {
        if (data == null) {
            data = this;
        } else if (data != this) {
            data.E(this);
        }
        return data;
    }

    /**
     * Returns a string formed from the dataInfo label and the array values.
     */
    public String toString() {
        return dataInfo.getLabel() + " " + x.toString();
    }
    
    private double[] x = new double[0];
    
    public static DataFactory getFactory(int n) {
        Factory factory = new Factory();
        factory.dataLength = n;
        return factory;
    }

    private static class Factory implements DataFactory, Serializable {
        
        int dataLength = 0;
        
        public Data makeData(String label, Dimension dimension) {
            DataDoubleArray data = new DataDoubleArray(label, dimension);
            data.setLength(dataLength);
            return data;
        }
        
        public Class getDataClass() {
            return double[].class;
        }
        
        public DataFactory copy() {
            return DataDoubleArray.getFactory(dataLength);
        }
        
    }

}