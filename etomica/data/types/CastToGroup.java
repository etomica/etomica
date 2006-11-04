package etomica.data.types;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.data.DataTag;
import etomica.data.types.DataGroup.DataInfoGroup;

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

    /**
     * Sole constructor.
     */
    public CastToGroup() {
    }

    protected DataInfo processDataInfo(DataInfo inputDataInfo) {
        Class inputClass = inputDataInfo.getClass();
        dataGroup = null;
        if (inputClass == DataGroup.class) {
            inputType = 0;
            dataGroup = null;
            return inputDataInfo;
        }
        inputType = 1;
        dataGroup = null;
        DataInfoGroup outputDataInfo = new DataInfoGroup(inputDataInfo.getLabel(), inputDataInfo.getDimension(), new DataInfo[]{inputDataInfo});
        return outputDataInfo;
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
            if (dataGroup == null) {
                dataGroup = new DataGroup(new Data[]{data});
            }
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
}
