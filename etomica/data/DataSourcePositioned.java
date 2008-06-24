package etomica.data;

import etomica.api.IBox;
import etomica.api.IVector;


/**
 * Interface for a DataSource that can return a value given an arbitrary
 * position.
 * 
 * @author Andrew Schultz
 */
public interface DataSourcePositioned {
    
    public void setBox(IBox box);
    
    public Data getData(IVector a);
    
    public IDataInfo getPositionDataInfo();
    
    public DataTag getTag();
}