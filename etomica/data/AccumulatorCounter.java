package etomica.data;

import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.units.Quantity;

/**
 * Data processor that simply counts the number of times its
 * <code>addData</code> method is invoked. Output is a DataDouble.
 */
public class AccumulatorCounter extends DataAccumulator {

    /**
     * @param parentElement
     * @param dataSource
     */
    public AccumulatorCounter() {
        dataInfo = new DataInfoDouble("Counter", Quantity.DIMENSION);
        data = new DataDouble();
    }
    
    /**
     * Returns null, indicating that any Data type is acceptable for input.
     */
    public DataProcessor getDataCaster(DataInfo incomingDataInfo) {
        return null;
    }

    /**
     * Does nothing.
     * 
     * @return the DataInfo for the output DataInteger
     */
    public DataInfo processDataInfo(DataInfo incomingDataInfo) {
        dataInfo.clearTags();
        dataInfo.addTags(incomingDataInfo.getTags());
        dataInfo.addTag(getTag());
        return dataInfo;
    }

    /**
     * Increments the counter. Argument is ignored.
     */
    protected void addData(Data dummyData) {
        data.x++;
    }

    /**
     * Returns the DataInteger with the count.
     */
    public Data getData() {
        return data;
    }

    /**
     * Sets count to zero.
     */
    public void reset() {
        data.x = 0;
    }

    private static final long serialVersionUID = 1L;

    protected final DataDouble data;
}
