/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.math.function.IFunction;


/**
 * Applies a simple scalar function to all elements of the Data
 * that passes through.  Function is set at construction and
 * cannot be subsequently changed.
 *
 * @author David Kofke
 */
public class DataProcessorFunction extends DataProcessor {

    public DataProcessorFunction(IFunction function) {
        this.function = function;
    }
    
    /**
     * Applies the function to all elements of the input data.
     */
    protected IData processData(IData inputData) {
        data.E(inputData);
        data.map(function);
        return data;
    }

    /**
     * Returns a copy of the given dataInfo
     */
    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        data = dataInfo.makeData();
        return dataInfo;
    }

    /**
     * Always returns null.
     */
    public DataPipe getDataCaster(IEtomicaDataInfo incomingDataInfo) {
        return null;
    }

    private static final long serialVersionUID = 1L;
    private final IFunction function;
    protected IData data;
}
