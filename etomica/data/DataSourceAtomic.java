package etomica.data;

import etomica.Atom;
import etomica.Data;
import etomica.DataSource;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */

/*
 * Created on Jul 23, 2004
 */

public interface DataSourceAtomic extends DataSource {
    
    public Data getData(Atom a);
    
}