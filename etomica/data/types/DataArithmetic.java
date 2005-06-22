package etomica.data.types;

import etomica.Data;
import etomica.utility.Function;

/**
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Style - Code Templates
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jun 15, 2005 by kofke
 */
public interface DataArithmetic extends DataNumeric {

    public void E(Data y);

    public void PE(DataArithmetic y);

    public void ME(DataArithmetic y);

    public void TE(DataArithmetic y);

    public void DE(DataArithmetic y);

    public void E(double y);

    public void PE(double y);

    public void TE(double y);

    public void map(Function function);

    public boolean isNaN();
}