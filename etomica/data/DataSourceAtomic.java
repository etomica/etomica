package etomica.data;

import etomica.atom.IAtom;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */

/*
 * Created on Jul 23, 2004
 */

public interface DataSourceAtomic extends DataSource {
    
    public Data getData(IAtom a);
    
    public IDataInfo getAtomDataInfo();
    
}