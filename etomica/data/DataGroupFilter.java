package etomica.data;

import etomica.data.types.CastToGroup;
import etomica.data.types.DataGroup;
import etomica.utility.Arrays;


/**
 * Extracts one or more of the Data instances wrapped in a DataGroup.  If a single instance,
 * output is the instance itself; if multiple instances, a new DataGroup with the subset is
 * output.
 *
 * @author David Kofke
 *
 */

public class DataGroupFilter extends DataProcessor {

    public DataGroupFilter(int index) {
        this(new int[] {index});
    }
    
    public DataGroupFilter(int[] indexes) {
        this.indexes = (int[])indexes.clone();
        singleInstance = (indexes.length == 1); 
    }

    public Data processData(Data data) {
        return outputData;
    }
    
    
    /* (non-Javadoc)
     * @see etomica.data.DataProcessor#makeOutputDataInfo(etomica.DataInfo)
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        Data[] inputData = ((DataGroup.Factory)inputDataInfo.getDataFactory()).getData();
        if(singleInstance) {
            if(indexes[0] < inputData.length) {
                outputData = inputData[indexes[0]];
            } else {
                throw new ArrayIndexOutOfBoundsException("DataFilter was constructed to extract a Data element with an index that is larger than the number of Data elements wrapped in the DataGroup. Number of elements: "+inputData.length+"; index array: "+Arrays.toString(indexes));
            }
        } else {
            Data[] pushedData = new Data[indexes.length];
            try {
                for (int i=0; i<indexes.length; i++) {
                    pushedData[i] = inputData[indexes[i]];
                }
            } catch(ArrayIndexOutOfBoundsException ex) {
                throw new ArrayIndexOutOfBoundsException("DataFilter was constructed to extract a Data element with an index that is larger than the number of Data elements wrapped in the DataGroup. Number of elements: "+inputData.length+"; index array: "+Arrays.toString(indexes));
            }
            outputData = new DataGroup(inputDataInfo.getLabel(),
                                            inputDataInfo.getDimension(),
                                            pushedData);
        }
        return outputData.getDataInfo();
    }

    /* (non-Javadoc)
     * @see etomica.DataSink#getDataCaster(etomica.DataInfo)
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        if(dataInfo.getDataClass() == DataGroup.class) {
            return null;
        }
        return new CastToGroup();
    }
    
    private Data outputData;
    private final boolean singleInstance;
    private final int[] indexes;
}
