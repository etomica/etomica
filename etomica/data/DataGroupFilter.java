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

    /**
     * Creates a filter that will take the DataGroup element indicated
     * by the index, counting from 0.
     */
    public DataGroupFilter(int index) {
        this(new int[] {index});
    }

    /**
     * Creates a filter that will take the set of elements from the DataGroup
     * corresponding to the given index values (with indexes counting from 0).
     */
    public DataGroupFilter(int[] indexes) {
        this.indexes = (int[])indexes.clone();
        singleInstance = (indexes.length == 1); 
    }

    /**
     * Returns the output data that was established with a previous call to 
     * processDataInfo.
     */
    public Data processData(Data data) {
        return outputData;
    }
    
    
    /**
     * Determines data that will be output, using the given DataInfo and
     * the index specification given at construction.  Returns the DataInfo
     * of the data that will be output.
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

    /**
     * Returns null if the given DataInfo is for a DataGroup; otherwise
     * returns a CastToGroup instance.
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
