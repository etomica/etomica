package etomica.data;

import etomica.api.IAtom;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */
public interface DataSourceAtomic {
    
    public Data getData(IAtom a);
    
    public IDataInfo getAtomDataInfo();
    
    public DataTag getTag();
}