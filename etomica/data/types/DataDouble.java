package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
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
public class DataDouble extends Data implements DataArithmetic {

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
     * Returns a copy of this instance.  Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataDouble(this);
    }

    public void E(Data y) {
        x = ((DataDouble)y).x;
    }

    public void E(double y) {
        x = y;
    }

    public void ME(DataArithmetic y) {
        x -= ((DataDouble)y).x;
    }

    public void PE(DataArithmetic y) {
        x += ((DataDouble)y).x;
    }

    public void TE(DataArithmetic y) {
        x *= ((DataDouble)y).x;
    }

    public void DE(DataArithmetic y) {
        x /= ((DataDouble)y).x;
    }

    public void PE(double y) {
        x += y;
    }

    public void TE(double y) {
        x *= y;
    }
    
    public boolean isNaN() {
        return Double.isNaN(x);
    }

    public void map(Function function) {
        x = function.f(x);
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
        return dataInfo.getLabel() + " " + Double.toString(x);
    }
    public double x;
    
    public interface Source extends DataSource {
        public DataDouble getDataDouble();
    }
}
