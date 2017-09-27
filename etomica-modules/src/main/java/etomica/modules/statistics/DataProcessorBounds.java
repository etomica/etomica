/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.statistics;

import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.units.dimensions.Null;

public class DataProcessorBounds extends DataProcessor {

    protected final DataDouble[] x;
    protected final DataDouble.DataInfoDouble[] xInfo;
    protected IDataSink[] sinks;
    protected DataTag[] tags;

    public DataProcessorBounds(IDataSink[] sinks) {
        super();
        x = new DataDouble[3];
        xInfo = new DataDouble.DataInfoDouble[3];
        String[] labels = new String[]{"xMinus", "x", "xPlus"};
        tags = new DataTag[3];
        for (int i = 0; i < 3; i++) {
            x[i] = new DataDouble();
            xInfo[i] = new DataDouble.DataInfoDouble(labels[i], Null.DIMENSION);
            tags[i] = new DataTag();
            xInfo[i].addTag(tags[i]);
            sinks[i].putDataInfo(xInfo[i]);
        }
        this.sinks = sinks;
    }

    public DataTag getTag(int i) {
        return tags[i];
    }

    @Override
    protected IData processData(IData inputData) {
        double avg = inputData.getValue(0);
        double err = inputData.getValue(1);
        for (int i = 0; i < 3; i++) {
            x[i].x = avg + (i - 1) * err;
            sinks[i].putData(x[i]);
        }

        return null;
    }

    @Override
    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        return null;
    }

    @Override
    public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
        return null;
    }
}
