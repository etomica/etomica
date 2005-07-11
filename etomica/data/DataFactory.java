package etomica.data;

import etomica.Data;
import etomica.units.Dimension;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jul 9, 2005 by kofke
 */
public interface DataFactory {

    public Data makeData(String label, Dimension dimension);
    
    public Class getDataClass();
    
    public DataFactory copy();

}
