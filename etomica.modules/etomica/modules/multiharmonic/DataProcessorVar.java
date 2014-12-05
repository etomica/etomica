/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.data.AccumulatorAverageCollapsingLog;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;

/**
 * Data Processor whose sole purpose in life is to discard the incoming data
 * and push on the standard deviation from the accumulator.  Using this means
 * that the standard deviation will be flowing down the data pipe at the same
 * time as the free energy to which it corresponds.
 *
 * @author Andrew Schultz
 */
public class DataProcessorVar extends DataProcessor {

    public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
        return null;
    }

    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        dataVar = inputDataInfo.makeData();
        return dataInfo;
    }

    protected IData processData(IData inputData) {
        dataVar.E(ac.getStdevLog());
        dataVar.TE(dataVar);
        return dataVar;
    }

    public void setAccumulator(AccumulatorAverageCollapsingLog newAc) {
        ac = newAc;
    }
    
    private static final long serialVersionUID = 1L;
    protected AccumulatorAverageCollapsingLog ac;
    protected IData dataVar;
}