/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import etomica.data.*;
import etomica.math.function.IFunction;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.dimensions.Dimension;

/**
 * Data object wrapping a single mutable value of type (Space) Tensor. Value is
 * final public and can be accessed directly. 
 * <p>
 * All arithmetic methods throw ClassCastException if given a Data instance that
 * is not of this type.
 * 
 * @author Andrew Schultz
 *  
 */

public class DataTensor implements IData {

    /**
     * Constructs a new instance with the given DataInfo, wrapping a new Tensor
     * instance from the given space.
     * 
     * @param space
     *            used to construct the wrapped Vector
     */
    public DataTensor(Space space) {
        super();
        x = space.makeTensor();
    }

    /**
     * Copies the elements of the given tensor (wrapped in the Data object) to
     * this vector.
     */
    public void E(IData y) {
        x.E(((DataTensor) y).x);
    }

    /**
     * Sets all tensor elements to the given value.
     */
    public void E(double y) {
        x.E(y);
    }

    /**
     * Minus-equals (-=) operation. Performed element-by-element.
     */
    public void ME(IData y) {
        x.ME(((DataTensor) y).x);
    }

    /**
     * Plus-equals (+=) operation. Performed element-by-element.
     */
    public void PE(IData y) {
        x.PE(((DataTensor) y).x);
    }

    /**
     * Times-equals (*=) operation. Performed element-by-element.
     */
    public void TE(IData y) {
        DataTensor yt = (DataTensor) y;
        for (int i = 0; i < x.D(); i++) {
            for (int j = 0; j < x.D(); j++) {
                x.setComponent(i, j, x.component(i, j) * yt.x.component(i, j));
            }
        }
    }

    /**
     * Divide-equals (/=) operation. Performed element-by-element.
     */
    public void DE(IData y) {
        x.DE(((DataTensor) y).x);
    }

    /**
     * Plus-equals (+=) operation. Adds given value to all elements.
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
     * Returns true if any vector element is not-a-number, as given by
     * Double.isNaN.
     */
    public boolean isNaN() {
        return x.isNaN();
    }

    /**
     * Maps the function on all vector elements, replacing each with the value
     * given by the function applied to it.
     */
    public void map(IFunction function) {
        x.map(function);
    }

    /**
     * Returns the number of elements in the wrapped tensor.
     */
    public int getLength() {
        return x.D() * x.D();
    }

    /**
     * Returns the i-th tensor value, indexing across the first row, then the
     * second, etc. For example, for a 2-D tensor, the indexes correspond to
     * 0=xx, 1=xy, 2=yx, 3=yy.
     */
    public double getValue(int i) {
        int D = x.D();
        if (i < 0 || i >= D*D)
            throw new IllegalArgumentException("Index out of bounds: " + i);
        int j = i / D;
        i -= D * j;
        return x.component(j, i);
    }

    /**
     * Returns a new array formed by the elements of the wrapped tensor.
     */
   public void assignTo(double[] array) {
        x.assignTo(array);
    }

    /**
     * Returns a string formed from the dataInfo label and the tensor values.
     */
    public String toString() {
        return x.toString();
    }

    /**
     * The wrapped tensor data.
     */
    public final Tensor x;
    
    public static class DataInfoTensor extends DataInfo {
        
        public DataInfoTensor(String label, Dimension dimension, Space space) {
            super(label, dimension);
            this.space = space;
        }
        
        public int getLength() {
            return space.D()*space.D();
        }
        
        public IDataInfoFactory getFactory() {
            return new DataInfoTensorFactory(this);
        }
        
        public Space getSpace() {
            return space;
        }
        
        public IData makeData() {
            return new DataTensor(space);
        }
        
        protected final Space space;
    }
    
    public static class DataInfoTensorFactory extends DataInfoFactory {
        protected DataInfoTensorFactory(DataInfoTensor template) {
            super(template);
            space = template.space;
        }
        
        public IDataInfo makeDataInfo() {
            DataInfoTensor dataInfo = new DataInfoTensor(label, dimension, space);
            DataTag[] tagArray = new DataTag[tags.size()];
            dataInfo.addTags((DataTag[])tags.toArray(tagArray));
            return dataInfo;
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
        
        protected Space space;
    }
}
