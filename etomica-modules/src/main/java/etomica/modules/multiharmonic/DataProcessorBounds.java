/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.data.*;

/**
 * Calculates upper or lower bound of the incoming data based on the standard
 * deviation obtained from the accumulator.
 *
 * @author Andrew Schultz
 */
public class DataProcessorBounds extends DataProcessor {

    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
        return null;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        data = inputDataInfo.makeData();
        return dataInfo;
    }

    protected IData processData(IData inputData) {
        if (inputData.getLength() < 2) {
            return data;
        }
        data.E(inputData);
        IData stdev = ac.getStdevLog();
        if (isUpperBound) {
            data.PE(stdev);
        }
        else {
            data.ME(stdev);
        }
        return data;
    }

    public void setAccumulator(AccumulatorAverageCollapsingLog newAc) {
        ac = newAc;
    }
    
    public void setIsUpperBound(boolean newIsUpperBound) {
        isUpperBound = newIsUpperBound;
    }
    
    private static final long serialVersionUID = 1L;
    protected IData data;
    protected IDataInfo dataInfo;
    protected AccumulatorAverageCollapsingLog ac;
    protected boolean isUpperBound;
}
