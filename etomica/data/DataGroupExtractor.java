package etomica.data;

import etomica.data.types.CastToGroup;
import etomica.data.types.DataGroup;


/**
 * Extracts a single Data instance from those wrapped in a DataGroup. 
 * Applies criterial definable by a DataJudge instance to assign a rank
 * to all elements in group, and returns the element having the lowest score.
 * If two or more elements have the minimum score, the first one encountered is the one returned.
 * Process applies recursively to DataGroups held by the DataGroup assigned for extraction. 
 *
 * @author David Kofke
 *
 */

public class DataGroupExtractor extends DataProcessor {

    public DataGroupExtractor(DataJudge judge) {
        this.judge = judge;
    }

    public Data processData(Data data) {
        return outputData;
    }
    
    
    /**
     * Find the most appropriate Data in the DataGroup.  Object with smallest (nearest zero)
     * rank assigned by DataJudge is selected.
     */
    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        MyWrapper outputWrapper = findBestData(inputDataInfo, 0);
        outputData = outputWrapper.data;
        return outputData.getDataInfo();
    }
    
    private MyWrapper findBestData(DataInfo inputDataInfo, int depth) {
        Data[] inputData = ((DataGroup.Factory)inputDataInfo.getDataFactory()).getData();
        MyWrapper outputWrapper = new MyWrapper(null, DataJudge.MAX_RANK);
        for(int i=0; i<inputData.length; i++) {
            Data testData = inputData[i];
            if(judge.doRecurse() && testData.getDataInfo().getClass() == DataGroup.class) {
                MyWrapper testWrapper = findBestData(testData.getDataInfo(), depth+1);
                if(testWrapper.rank < outputWrapper.rank) {
                    outputWrapper.data = testWrapper.data;
                    outputWrapper.rank = testWrapper.rank;
                }
            } else {
                int rank = judge.rankData(testData, i, depth);
                if(rank < outputWrapper.rank) {
                    outputWrapper.data = testData;
                    outputWrapper.rank = rank;
                }
            }
        }
        return outputWrapper;
    }

    /* (non-Javadoc)
     * @see etomica.DataSink#getDataCaster(etomica.DataInfo)
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        if(dataInfo.getClass() == DataGroup.class) {
            return null;
        }
        return new CastToGroup();
    }
    
    private Data outputData;
    private final DataJudge judge;
 
    //simple data structure holding Data and its rank
    private static class MyWrapper {
        Data data;
        int rank;
        MyWrapper(Data data, int rank) {
            this.data = data;
            this.rank = rank;
        }
    }
}