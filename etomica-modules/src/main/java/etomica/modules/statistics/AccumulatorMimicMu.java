/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;

import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorBox;

public class AccumulatorMimicMu extends AccumulatorAverageFixed {

    protected final IntegratorBox integrator;
    protected DataGroup data;

    public AccumulatorMimicMu(IntegratorBox integrator) {
        super();
        this.integrator = integrator;
    }

    public void putDataInfo(IDataInfo dataInfo) {
        dataSink.putDataInfo(dataInfo);
        data = (DataGroup) dataInfo.makeData();
    }

    public IData processData(IData inputData) {
        double temperature = integrator.getTemperature();
        DataGroup g = (DataGroup) inputData;
        DataDouble avg = (DataDouble) g.getData(AVERAGE.index);
        DataDouble err = (DataDouble) g.getData(ERROR.index);
        DataDouble stdev = (DataDouble) g.getData(STANDARD_DEVIATION.index);
        DataDouble avgOut = (DataDouble) data.getData(AVERAGE.index);
        DataDouble errOut = (DataDouble) data.getData(ERROR.index);
        DataDouble stdevOut = (DataDouble) data.getData(STANDARD_DEVIATION.index);
        errOut.x = temperature * err.x / avg.x;
        stdevOut.x = temperature * stdev.x / avg.x;
        avgOut.x = -temperature * Math.log(avg.x);
        return data;
    }
}
