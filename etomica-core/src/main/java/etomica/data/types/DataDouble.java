/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import etomica.math.function.IFunction;
import etomica.data.DataInfo;
import etomica.data.DataInfoFactory;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataInfoFactory;
import etomica.units.Dimension;

/**
 * Data object wrapping a single mutable value of type <tt>double</tt>. Value is
 * public and can be accessed directly.
 * <p>
 * All arithmetic methods throw ClassCastException if given a Data instance that
 * is not of this type.
 * 
 * @author David Kofke and Andrew Schultz
 *  
 */
public class DataDouble implements IData, java.io.Serializable {

    /**
     * Constructs a new instance with given descriptors.
     */
    public DataDouble() {
        super();
    }

    /**
     * Sets the wrapped double to the value in the given instance.
     */
    public void E(IData y) {
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
    public void ME(IData y) {
        x -= ((DataDouble) y).x;
    }

    /**
     * Plus-equals (+=) operation. Adds the value in the given instance to this
     * instance's value.
     */
    public void PE(IData y) {
        x += ((DataDouble) y).x;
    }

    /**
     * Times-equals (*=) operation. Replaces the value in this instance with its
     * value times the value in the given instance.
     */
    public void TE(IData y) {
        x *= ((DataDouble) y).x;
    }

    /**
     * Divide-equals (/=) operation. Divides this value by the value in the
     * given instance.
     */
    public void DE(IData y) {
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
     * Times-equals (*=) operation. Multiplies this value by that in the given
     * instance.
     */
    public void TE(double y) {
        x *= y;
    }

    /**
     * Divide-equals (/=) operation. Divides this value by that in the given
     * instance.
     */
    public void DE(double y) {
        x /= y;
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
    public void map(IFunction function) {
        x = function.f(x);
    }

    /**
     * Returns a new one-element array formed from the current value of this
     * instance.
     */
    public void assignTo(double[] array) {
        array[0] = x;
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

    
    /**
     * Returns a string formed from the dataInfo label and this value.
     */
    public String toString() {
        return Double.toString(x);
    }

    private static final long serialVersionUID = 1L;
    /**
     * The wrapped data value held by this object.
     */
    public double x;
    
    public static class DataInfoDouble extends DataInfo {
        public DataInfoDouble(String label, Dimension dimension) {
            super(label, dimension);
        }
        
        public int getLength() {
            return 1;
        }
        
        public IEtomicaDataInfoFactory getFactory() {
            return new DataInfoDoubleFactory(this);
        }
        
        public IData makeData() {
            return new DataDouble();
        }

        private static final long serialVersionUID = 1L;
    }
    
    public static class DataInfoDoubleFactory extends DataInfoFactory {
        protected DataInfoDoubleFactory(DataInfoDouble template) {
            super(template);
        }
        
        public IEtomicaDataInfo makeDataInfo() {
            DataInfoDouble dataInfo = new DataInfoDouble(label, dimension);
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
        }

        private static final long serialVersionUID = 1L;
    }
}
