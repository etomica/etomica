package etomica.data;

import etomica.api.IAtom;
import etomica.api.IData;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */
public interface DataSourceAtomic {
    
    public IData getData(IAtom a);
    
    public IEtomicaDataInfo getAtomDataInfo();
    
    public DataTag getTag();
}