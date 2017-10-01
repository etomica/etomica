/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * 
 */
package etomica.modules.multiharmonic;

import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataSourceScalar;
import etomica.data.IData;
import etomica.data.IDataInfo;

/**
 * Data Processor which calculates the bias in the incoming data based on the
 * true value of the free energy calculated from the 
 *
 * @author Andrew Schultz
 */
public class DataProcessorBias extends DataProcessor {

    public DataProcessorBias(DataSourceScalar trueDataSource) {
        this.trueDataSource = trueDataSource;
    }
    
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
        if (inputData.getLength() < 1) {
            return data;
        }

        data.E(inputData);
        data.PE(-trueDataSource.getDataAsScalar());
        return data;
    }

    private static final long serialVersionUID = 1L;
    protected IData data;
    protected IDataInfo dataInfo;
    protected DataSourceScalar trueDataSource;
}
