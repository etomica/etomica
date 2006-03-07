package etomica.modules.multiharmonic;


import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.data.types.DataArray;
import etomica.units.Dimension;

/**
 * A DataProcessor that unwraps a DataArray, returning the first element.
 * 
 * @author Andrew Schultz
 */
public class CastArrayUnwrap extends DataProcessor {

    /**
     * Sole constructor.
     */
    public CastArrayUnwrap() {
    }
    
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        if (inputDataInfo.getDataClass() != DataArray.class) {
            throw new IllegalArgumentException("can cast only from DataArray");
        }
        String label = inputDataInfo.getLabel();
        Dimension dimension = inputDataInfo.getDimension();
        DataArray.Factory factory = (DataArray.Factory)inputDataInfo.getDataFactory();
        DataFactory subFactory = factory.getArrayElementFactory();
        
        return subFactory.makeData(label,dimension).getDataInfo();
    }
    
    /**
     * Returns the first element from the given DataArray
     */
    protected Data processData(Data data) {
        DataArray dataArray = (DataArray)data;
        return dataArray.getData(0);
    }
    
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }
}
