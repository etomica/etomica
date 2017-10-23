/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.integrator.IntegratorBox;

/**
 * DataProcessor that takes in <exp(-beta*deltaU)> and spits out
 * mu = -T*ln<exp(-beta*deltaU)>
 */
public class DataProcessorMu extends DataProcessor {

    protected final IntegratorBox integrator;
    protected final DataDouble data;

    public DataProcessorMu(IntegratorBox integrator) {
        super();
        this.integrator = integrator;
        data = new DataDouble();
    }

    public IDataInfo processDataInfo(IDataInfo dataInfo) {
        if (dataSink != null) dataSink.putDataInfo(dataInfo);
        return dataInfo;
    }

    public IData processData(IData inputData) {
        double temperature = integrator.getTemperature();
        double x = inputData.getValue(0);
        data.x = -temperature * Math.log(x);
        return data;
    }
}
