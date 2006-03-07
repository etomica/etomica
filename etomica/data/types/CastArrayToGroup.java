package etomica.data.types;


import etomica.data.Data;
import etomica.data.DataFactory;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.units.Dimension;

/**
 * A DataProcessor that converts a DataArray instance into a DataGroup.
 * If the DataArray is multidimensional, it will be flattened into a 1-D array.
 *
 * @author Andrew Schultz
 */
public class CastArrayToGroup extends DataProcessor {

    /**
     * Sole constructor.
     */
    public CastArrayToGroup() {
    }
    
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        if (inputDataInfo.getDataClass() != DataArray.class) {
            throw new IllegalArgumentException("can cast only from DataArray");
        }
        String label = inputDataInfo.getLabel();
        Dimension dimension = inputDataInfo.getDimension();
        DataArray.Factory factory = (DataArray.Factory)inputDataInfo.getDataFactory();
        DataFactory subFactory = factory.getArrayElementFactory();
        
        Data[] groupData = new Data[factory.getArrayLength()];
        for (int i=0; i<groupData.length; i++) {
            groupData[i] = subFactory.makeData(label, dimension);
        }
        // this is bogus, but the DataInfo is ok.  We'll fix it later
        outputData = new DataGroup(label, groupData);

        return outputData.getDataInfo();
    }
    
    /**
     * Extracts a double from the input data and returns it encapsulated in a
     * DataDouble.
     * 
     * @return a DataGroup holding the Data elements from the given DataArray
     */
    protected Data processData(Data data) {
        DataArray dataArray = (DataArray)data;
        boolean needUpdate = false;
        if (dataArray.getLength() != outputData.getNData()) {
            needUpdate = true;
        }
        else {
            for (int i=0; i<dataArray.getLength(); i++) {
                if (dataArray.getData(i) != outputData.getData(i)) {
                    needUpdate = true;
                    break;
                }
            }
        }
        if (needUpdate) {
            // something has changed.  We need to rewrap the array of Data in 
            // a new DataGroup
            Data[] groupData = new Data[dataArray.getLength()];
            for (int i=0; i<groupData.length; i++) {
                groupData[i] = dataArray.getData(i);
            }
            outputData = new DataGroup(dataArray.getDataInfo().getLabel(),groupData);
        }
        return outputData;
    }
    
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }

    private DataGroup outputData;
}
