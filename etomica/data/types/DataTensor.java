package etomica.data.types;

import java.io.Serializable;

import etomica.Space;
import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.space.Tensor;
import etomica.units.Dimension;
import etomica.utility.Function;

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

/*
 * History Created on Jun 15, 2005 by andrew
 */
public class DataTensor extends Data implements DataArithmetic {

    /**
     * Constructs a new instance with the given DataInfo, wrapping a new Vector
     * instance from the given space.
     * 
     * @param space
     *            used to construct the wrapped Vector
     * @param label
     *            a descriptive label for this data
     * @param dimension
     *            the physical dimensions of the data
     */
    public DataTensor(Space space, String label, Dimension dimension) {
        super(new DataInfo(label, dimension, getFactory(space)));
        x = space.makeTensor();
    }

    /**
     * Copy constructor.
     */
    public DataTensor(DataTensor data) {
        super(data);
        x = (Tensor) data.x.clone();
    }

    /**
     * Returns a deep copy of this instance. Returned object has its own
     * instances of all fields, set equal to the values of this instance's
     * fields. This instance and the copy share the same DataInfo, which is immutable.
     */
    public Data makeCopy() {
        return new DataTensor(this);
    }

    /**
     * Copies the elements of the given tensor (wrapped in the Data object) to
     * this vector.
     */
    public void E(Data y) {
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
    public void ME(DataArithmetic y) {
        x.ME(((DataTensor) y).x);
    }

    /**
     * Plus-equals (+=) operation. Performed element-by-element.
     */
    public void PE(DataArithmetic y) {
        x.PE(((DataTensor) y).x);
    }

    /**
     * Times-equals (*=) operation. Performed element-by-element.
     */
    public void TE(DataArithmetic y) {
        x.TE(((DataTensor) y).x);
    }

    /**
     * Divide-equals (/=) operation. Performed element-by-element.
     */
    public void DE(DataArithmetic y) {
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
    public void map(Function function) {
        x.map(function);
    }

    /**
     * Returns the number of elements in the wrapped vector.
     */
    public int getLength() {
        return x.length();
    }

    /**
     * Returns the i-th tensor value, indexing across the first row, then the
     * second, etc. For example, for a 2-D tensor, the indexes correspond to
     * 0=xx, 1=xy, 2=yx, 3=yy.
     */
    public double getValue(int i) {
        int D = x.length();
        if (i < 0 || i >= D)
            throw new IllegalArgumentException("Index out of bounds: " + i);
        int j = i / D;
        i -= D * j;
        return x.component(j, i);
    }

    /**
     * Returns a new array formed by the elements of the wrapped vector.
     */
   public void assignTo(double[] array) {
        x.assignTo(array);
    }

    /**
     * Returns a string formed from the dataInfo label and the tensor values.
     */
    public String toString() {
        return dataInfo.getLabel() + " " + x.toString();
    }
    
    /**
     * Returns a DataFactory that makes a DataTensor for the given Space.
     */
    public static DataFactory getFactory(Space space) {
        if(FACTORY == null || FACTORY.space != space) { 
            FACTORY = new Factory(space);
        }
        return FACTORY;
    }

    /**
     * The wrapped tensor data.
     */
    public final Tensor x;
    
    private transient static Factory FACTORY = null;
    
    /**
     * DataFactory that constructs DataTensor instances of a given shape.
     * Instantiate using the static DataTensor.getFactory method. 
     */
    public static class Factory implements DataFactory, Serializable {
        
        private final Space space;
        
        Factory(Space space) {
            this.space = space;
        }
        
        /**
         * Makes and returns a new DataTensor using the given label and dimension
         * for its DataInfo.
         */
        public Data makeData(String label, Dimension dimension) {
            return new DataTensor(space, label, dimension);
        }
        
        /**
         * Returns the class of the data made by this factory.
         * @return DataTensor.class
         */
        public Class getDataClass() {
            return DataTensor.class;
        }
        
        /**
         * @return the space making the wrapped tensor
         */
        public Space getSpace() {
            return space;
        }
    }//end of Factory

}
