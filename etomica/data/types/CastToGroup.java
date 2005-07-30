package etomica.data.types;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.DataProcessor;

/**
 * A DataProcessor that wraps a Data instance into a DataGroup.  Has no effect
 * if input Data is already a DataGroup. 
 *
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 21, 2005 by kofke
 */
public class CastToGroup extends DataProcessor {

    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        Class inputClass = inputDataInfo.getClass();
        dataGroup = null;
        if (inputClass == DataGroup.class) {
            inputType = 0;
            return inputDataInfo;
        } 
        //wrapping; it is necessary to make the DataGroup now, becuase we need it to return the DataInfo.
        //its DataInfo cannot be constructed without the wrapped Data, so we need to wrap a new instance
        //of the Data and then copy the data to it during processData
        inputType = 1;
        wrappedData = inputDataInfo.getDataFactory().makeData(inputDataInfo.getLabel(), inputDataInfo.getDimension());
        dataGroup = new DataGroup(inputDataInfo.getLabel()+" wrapped", inputDataInfo.getDimension(),
                        new Data[] {wrappedData});
        return dataGroup.getDataInfo();
    }
    
    /**
     * Extracts a double from the input data and returns it encapsulated in a
     * DataDouble.
     * 
     * @param data
     *            a Data instance of the type indicated by the DataInfo at
     *            construction
     * @return a DataDouble holding the value cast from the given Data; the same
     *         instance is returned with every invocation.
     * 
     * @throws ClassCastException
     *             if the given Data is not of the same type as indicated by the
     *             DataInfo given at construction
     */
    protected Data processData(Data data) {
        switch (inputType) {
        case 0:
            return data;
        case 1:
            wrappedData.E(data);
            return dataGroup;
        default:
            throw new Error("Assertion error.  Input type out of range: "+inputType);
        }
    }
    
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }

    private DataGroup dataGroup;
    private int inputType;
    private Data wrappedData;
}
