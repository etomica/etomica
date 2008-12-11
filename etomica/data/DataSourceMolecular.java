package etomica.data;

import etomica.api.IData;
import etomica.api.IMolecule;


/**
 * Interface for a DataSource that can return a value given an arbitrary atom.
 */
public interface DataSourceMolecular {
    
    public IData getData(IMolecule a);
    
    public IEtomicaDataInfo getMoleculeDataInfo();
    
    public DataTag getTag();
}