package etomica.data;

import etomica.api.IAtomLeaf;
import etomica.api.IData;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */
public interface DataSourceAtomic {
    
    public IData getData(IAtomLeaf a);
    
    public IEtomicaDataInfo getAtomDataInfo();
    
    public DataTag getTag();
}