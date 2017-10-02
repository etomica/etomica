/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.types;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
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
public class CastToGroup extends DataProcessor {

    private DataGroup dataGroup;
    private int inputType;
    
    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        Class inputClass = inputDataInfo.getClass();
        dataGroup = null;
        if (inputClass == DataGroup.class) {
            inputType = 0;
            dataGroup = null;
            return inputDataInfo;
        }
        inputType = 1;
        dataGroup = null;
        DataInfoGroup outputDataInfo = new DataInfoGroup(inputDataInfo.getLabel(), inputDataInfo.getDimension(), new IDataInfo[]{inputDataInfo});
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
    protected IData processData(IData data) {
        switch (inputType) {
        case 0:
            return data;
        case 1:
            if (dataGroup == null) {
                dataGroup = new DataGroup(new IData[]{data});
            }
            return dataGroup;
        default:
            throw new Error("Assertion error.  Input type out of range: "+inputType);
        }
    }
}
