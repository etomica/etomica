package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;

/**
 * A DataProcessor that effectively wraps a Data instance into a DataGroup.  Has no effect
 * if input Data is already a DataGroup.
 * <p>
 * A new instance of the input Data is wrapped in the output DataGroup, and the
 * processData method copies the input values to those in the copy.
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
        dataGroup = new DataGroup(inputDataInfo.getLabel()+" wrapped", new Data[] {wrappedData});
        return dataGroup.getDataInfo();
    }
    
    /**
     * Processes the input Data to update the output DataGroup.  If the input is
     * a DataGroup, is is simply returned; otherwise it values are copied to the
     * wrapped Data, and the wrapping DataGroup is returned.
     * 
     * @param data
     *            a Data instance of the type indicated by the DataInfo at
     *            the most recent call to processDataInfo
     * @return a DataGroup holding the values from given Data
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
    
    /**
     * Returns null, indicating that this DataProcessor can accept any Data type.
     */
    public DataProcessor getDataCaster(DataInfo info) {
        return null;
    }

    private DataGroup dataGroup;
    private int inputType;
    private Data wrappedData;
}
