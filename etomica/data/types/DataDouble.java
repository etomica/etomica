package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.utility.Function;

/**
 * Data object wrapping a single mutable value of type double. Value is
 * public and can be accessed directly.<br>
 * All arithmetic methods throw ClassCastException if given a Data instance that
 * is not of this type.
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */

/*
 * History Created on Jun 15, 2005 by kofke
 */
public class DataDouble extends Data implements DataArithmetic {

    /**
     * Constructs a new instance with the given DataInfo.
     */
    public DataDouble(DataInfo dataInfo) {
        super(dataInfo);
    }

    /**
     * Copy constructor.
     */
    public DataDouble(DataDouble data) {
        super(data);
        x = data.x;
    }

    /**
     * Returns a copy of this instance. Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataDouble(this);
    }

    /**
     * Sets the wrapped double to the value in the given instance.
     */
    public void E(Data y) {
        x = ((DataDouble) y).x;
    }

    /**
     * Sets the wrapped double to the given value.
     */
    public void E(double y) {
        x = y;
    }

    /**
     * Minus-equals (-=) operation. Subtracts the value in the given instance
     * from this instance's value.
     */
    public void ME(DataArithmetic y) {
        x -= ((DataDouble) y).x;
    }

    /**
     * Plus-equals (+=) operation. Adds the value in the given instance to this
     * instance's value.
     */
    public void PE(DataArithmetic y) {
        x += ((DataDouble) y).x;
    }

    /**
     * Times-equals (*=) operation. Replaces the value in this instance with its
     * value times the value in the given instance.
     */
    public void TE(DataArithmetic y) {
        x *= ((DataDouble) y).x;
    }

    /**
     * Divide-equals (/=) operation. Divides this value by the value in the
     * given instance.
     */
    public void DE(DataArithmetic y) {
        x /= ((DataDouble) y).x;
    }

    /**
     * Plus-equals (+=) operation. Adds the given value to this instance's
     * value.
     */
    public void PE(double y) {
        x += y;
    }

    /**
     * Plus-equals (*=) operation. Multiplies this value by that in the given
     * instance.
     */
    public void TE(double y) {
        x *= y;
    }

    /**
     * Returns true if the this instance's value is NaN.
     */
    public boolean isNaN() {
        return Double.isNaN(x);
    }

    /**
     * Maps the function on this instance's value, replacing it with the value
     * returned by the function.
     */
    public void map(Function function) {
        x = function.f(x);
    }

    /**
     * Returns a new one-element array formed from the current value of this
     * instance.
     */
    public double[] toArray() {
        return new double[] { x };
    }

    /**
     * Returns 1, indicating that this data object holds one value.
     */
    public int getLength() {
        return 1;
    }

    /**
     * Returns this instance's value if argument is 0, otherwise throw exception
     * 
     * @throws IllegalArgumentException
     *             if i != 0
     */
    public double getValue(int i) {
        if (i == 0) {
            return x;
        }
        throw new IllegalArgumentException(
                "Only permissible value for index is 0; value given is " + i);
    }

    
    public DataArithmetic toArithmetic(DataArithmetic data) {
        if (data == null) {
            data = this;
        } else if (data != this) {
            data.E(this);
        }
        return this;
    }

    /**
     * Returns a string formed from the dataInfo label and this value.
     */
    public String toString() {
        return dataInfo.getLabel() + " " + Double.toString(x);
    }

    /**
     * The wrapped data value held by this object.
     */
    public double x;

    /**
     * Interface indicating the a data source gives data of type DataDouble.
     */
    public interface Source extends DataSource {

        public DataDouble getDataDouble();
    }
}
