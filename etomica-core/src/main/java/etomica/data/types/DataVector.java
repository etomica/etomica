/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import etomica.data.*;
import etomica.math.function.IFunction;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;

/**
 * Data object wrapping a single mutable value of type (Space) Vector. Value is
 * final public and can be accessed directly. 
 * <p>
 * All arithmetic methods throw ClassCastException if given a Data instance that
 * is not of this type.
 * 
 * @author David Kofke
 *  
 */
public class DataVector implements IData, java.io.Serializable {

    /**
     * Constructs a new instance with the given DataInfo, wrapping a new Vector
     * instance from the given space.
     * 
     * @param space
     *            used to construct the wrapped Vector
     */
    public DataVector(Space space) {
        super();
        x = space.makeVector();
    }

    /**
     * Copies the elements of the given vector (wrapped in the Data object)
     * to this vector.
     */
    public void E(IData y) {
        x.E(((DataVector) y).x);
    }

    /**
     * Sets all vector elements to the given value.
     */
    public void E(double y) {
        x.E(y);
    }

    /**
     * Minus-equals (-=) operation.  Performed element-by-element.
     */
    public void ME(IData y) {
        x.ME(((DataVector) y).x);
    }

    /**
     * Plus-equals (+=) operation. Performed element-by-element.
     */
    public void PE(IData y) {
        x.PE(((DataVector) y).x);
    }

    /**
     * Times-equals (*=) operation. Performed element-by-element.
     */
    public void TE(IData y) {
        x.TE(((DataVector) y).x);
    }

    /**
     * Divide-equals (/=) operation. Performed element-by-element.
     */
    public void DE(IData y) {
        x.DE(((DataVector) y).x);
    }

    /**
     * Plus-equals (+=) operation.  Adds given value to all elements.
     */
    public void PE(double y) {
        x.PE(y);
    }

    /**
     * Times-equals (*=) operation. Multiplies all elements by the given value.
     */
    public void TE(double y) {
        x.TE(y);
    }

    /**
     * Divide-equals (*=) operation. Divides all elements by the given value.
     */
    public void DE(double y) {
        x.TE(1.0/y);
    }

    /**
     * Returns true if any vector element is not-a-number, as given by Double.isNaN.
     */
    public boolean isNaN() {
        return x.isNaN();
    }

    /**
     * Maps the function on all vector elements, replacing each with the
     * value given by the function applied to it.
     */
    public void map(IFunction function) {
        x.map(function);
    }

    /**
     * Returns the number of elements in the wrapped vector.
     */
    public int getLength() {
        return x.getD();
    }

    /**
     * Returns the i-th vector value.
     */
    public double getValue(int i) {
        if(i < 0 || i>= x.getD()) throw new IllegalArgumentException("Illegal value: " + i);
        return x.getX(i);
    }

    /**
     * Assigns the elements of the wrapped vector to the given array.
     */
    public void assignTo(double[] array) {
        x.assignTo(array);
    }

    /**
     * Returns a string formed from the dataInfo label and the vector values.
     */
    public String toString() {
        return x.toString();
    }
    
    private static final long serialVersionUID = 1L;
    /**
     * The wrapped vector data.
     */
    public final Vector x;
    
    public static class DataInfoVector extends DataInfo {
        
        public DataInfoVector(String label, Dimension dimension, Space space) {
            super(label, dimension);
            this.space = space;
        }
        
        public int getLength() {
            return space.D();
        }
        
        public IDataInfoFactory getFactory() {
            return new DataInfoVectorFactory(this);
        }
        
        public Space getSpace() {
            return space;
        }
        
        public IData makeData() {
            return new DataVector(space);
        }

        private static final long serialVersionUID = 1L;
        protected final Space space;
    }

    public static class DataInfoVectorFactory extends DataInfoFactory {
        protected DataInfoVectorFactory(DataInfoVector template) {
            super(template);
            space = template.space;
        }
        
        public IDataInfo makeDataInfo() {
            return new DataInfoVector(label, dimension, space);
        }
        
        /**
         * Sets the Space
         */
        public void setSpace(Space newSpace) {
            space = newSpace;
        }
        
        /**
         * Returns the Space
         */
        public Space getSpace() {
            return space;
        }

        private static final long serialVersionUID = 1L;
        protected Space space;
    }
}
