/*
 * Created on Jul 23, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.data;

import etomica.Atom;
import etomica.Data;
import etomica.DataSource;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */
public interface DataSourceAtomic extends DataSource {
    
    public Data getData(Atom a);
    
}