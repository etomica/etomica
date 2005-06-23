package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.Space;
import etomica.space.Tensor;
import etomica.utility.Function;

/*
 * History
 * Created on Jun 15, 2005 by andrew
 */
public class DataTensor extends Data implements DataArithmetic {

    public DataTensor(Space space, DataInfo dataInfo) {
        super(dataInfo);
        x = space.makeTensor();
    }

    /**
     * Copy constructor.
     */
    public DataTensor(DataTensor data) {
        super(data);
        x = (Tensor)data.x.clone();
    }
    
    /**
     * Returns a copy of this instance.  Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataTensor(this);
    }

    public void E(Data y) {
        x.E(((DataTensor)y).x);
    }

    public void E(double y) {
        x.E(y);
    }

    public void ME(DataArithmetic y) {
        x.ME(((DataTensor)y).x);
    }

    public void PE(DataArithmetic y) {
        x.PE(((DataTensor)y).x);
    }

    public void TE(DataArithmetic y) {
        x.TE(((DataTensor)y).x);
    }

    public void DE(DataArithmetic y) {
        x.DE(((DataTensor)y).x);
    }

    public void PE(double y) {
        x.PE(y);
    }

    public void TE(double y) {
        x.TE(y);
    }
    
    public boolean isNaN() {
        return x.isNaN();
    }

    public void map(Function function) {
        x.map(function);
    }
    
    public int getLength() {
        return x.length();
    }
    
    //XXX this doesn't work
    public double getValue(int i) {
        throw new RuntimeException("Method not appropriate for Tensor");
//        return x.x(i);
    }
    
    //XXX this doesn't work
    public double[] toArray() {
        throw new RuntimeException("Method not appropriate for Tensor");
//        return x.toArray();
    }
    
    public DataArithmetic toArithmetic(DataArithmetic data) {
        if (data == null) {
            data = this;
        }
        else if (data != this) {
            data.E(this);
        }
        return this;
    }
    
    public String toString() {
        return dataInfo.getLabel() + " " + x.toString();
    }
    public final Tensor x;
}
