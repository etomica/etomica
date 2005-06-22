package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.Space;
import etomica.space.Vector;
import etomica.utility.Function;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 15, 2005 by kofke
 */
public class DataVector extends Data implements DataArithmetic {

    public DataVector(Space space, DataInfo dataInfo) {
        super(dataInfo);
        x = space.makeVector();
    }

    /**
     * Copy constructor.
     */
    public DataVector(DataVector data) {
        super(data);
        x = (Vector)data.x.clone();
    }
    
    /**
     * Returns a copy of this instance.  Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataVector(this);
    }

    public void E(Data y) {
        x.E(((DataVector)y).x);
    }

    public void E(double y) {
        x.E(y);
    }

    public void ME(DataArithmetic y) {
        x.ME(((DataVector)y).x);
    }

    public void PE(DataArithmetic y) {
        x.PE(((DataVector)y).x);
    }

    public void TE(DataArithmetic y) {
        x.TE(((DataVector)y).x);
    }

    public void DE(DataArithmetic y) {
        x.DE(((DataVector)y).x);
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
    public final Vector x;
}
