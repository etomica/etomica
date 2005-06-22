package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;


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
public class DataInteger extends Data implements DataNumeric {

    public DataInteger(DataInfo dataInfo) {
        super(dataInfo);
    }

    /**
     * Copy constructor.
     */
    public DataInteger(DataInteger data) {
        super(data);
        x = data.x;
    }
    
    /**
     * Returns a copy of this instance.  Returned object has its own instances of
     * all fields, set equal to the values of this instance's fields.
     */
    public Data makeCopy() {
        return new DataInteger(this);
    }

    public void E(Data y) {
        x = ((DataInteger)y).x;
    }

    public void E(int y) {
        x = y;
    }
    
    public DataArithmetic toArithmetic(DataArithmetic data) {
        if (data == null) {
            data = new DataDouble(getDataInfo());
        }
        ((DataDouble)data).x = x;
        return data;
    }

    public String toString() {
        return dataInfo.getLabel() + " " + Integer.toString(x);
    }
    public int x;
}
