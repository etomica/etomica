package etomica.models.hexane;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.data.DataTag;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;

/**
 * Takes a N dimensional array and flattens it out by removing the first dimension 
 * and creating an N-1 dimensional array.
 * @author nancycribbin
 *
 */
public class DataProcessorArrayFlatten extends DataProcessor {

    /**
     * Constructor that doesn't do anything.
     */
    public DataProcessorArrayFlatten(){
        tag = new DataTag();
    }
    
    protected Data processData(Data inputData) {
        outputData = new DataDoubleArray(shapeNew, ((DataDoubleArray)inputData).getData());
        return outputData;
    }

    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        if(!(inputDataInfo instanceof DataInfoDoubleArray)){
            throw new IllegalArgumentException("DataProcessorArrayFlatten needs DataDoubleArray argument");
        }
        //Get the shape of the data we expect to come in, reduce it by one, then 
        // create a new DataDoubleArray of that shape.
        int[] shapeOld = ((DataInfoDoubleArray)inputDataInfo).getArrayShape();
        shapeNew = new int[shapeOld.length-1];
        shapeNew[0] = shapeOld[0] * shapeOld[1];
        for(int i = 1; i < shapeNew.length; i++){
            shapeNew[i] = shapeOld[i+1];
        }
        
        outputDataInfo = new DataInfoDoubleArray("unlabeled", inputDataInfo.getDimension(), shapeNew);
        
        return outputDataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * returns null (non-Javadoc)
     */
    public DataProcessor getDataCaster(DataInfo dataInfo) {
        return null;
    }

    private DataDoubleArray outputData;
    private DataInfoDoubleArray outputDataInfo;
    private int[] shapeNew;
    protected final DataTag tag;
}
