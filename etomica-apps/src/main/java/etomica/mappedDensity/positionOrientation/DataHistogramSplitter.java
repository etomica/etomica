/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.mappedDensity.positionOrientation;

import etomica.data.*;
import etomica.data.types.DataFunction;
import etomica.units.dimensions.Null;

public class DataHistogramSplitter implements IDataSink {

    protected DataTag[] tags;
    protected IDataSink[] sinks;
    protected DataFunction.DataInfoFunction[] dataInfo;
    protected DataFunction[] data;

    public DataHistogramSplitter() {
        tags = new DataTag[0];
        sinks = new IDataSink[0];
    }

    public void setDataSink(int i, IDataSink dataSink) {
        sinks[i] = dataSink;
        if (dataInfo[i]!=null) sinks[i].putDataInfo(dataInfo[i]);
    }

    public DataTag getTag(int i) {
        return tags[i];
    }

    @Override
    public void putData(IData inData) {
        for (int i=0; i<sinks.length; i++) {
            if (sinks[i] == null) continue;
            double[] y = data[i].getData();
            ((DataFunction)inData).assignColumnTo(i, y);
            sinks[i].putData(data[i]);
        }
    }

    @Override
    public void putDataInfo(IDataInfo incomingDataInfo) {
        DataFunction.DataInfoFunction inDataInfo = (DataFunction.DataInfoFunction)incomingDataInfo;
        int[] shape = inDataInfo.getArrayShape();
        if (shape[0] * shape[1] == 0) return;
        int nHistograms = shape[0];
        DataSourceIndependent inDsi = inDataInfo.getXDataSource();
        if (nHistograms == sinks.length) {
            if (shape[1] == dataInfo[0].getLength()) return;
            for (int i=0; i<sinks.length; i++) {
                dataInfo[i] = new DataFunction.DataInfoFunction("stuff", Null.DIMENSION, new DataSourceIndependentSimple(inDsi.getIndependentData(1).getData(),inDsi.getIndependentDataInfo(1)));
                dataInfo[i].addTag(tags[i]);
                data[i] = new DataFunction(new int[]{shape[1]});
                sinks[i].putDataInfo(dataInfo[i]);
            }
            return;
        }
        sinks = new IDataSink[nHistograms];
        tags = new DataTag[nHistograms];
        dataInfo = new DataFunction.DataInfoFunction[nHistograms];
        data = new DataFunction[nHistograms];
        for (int i=0; i<tags.length; i++) {
            tags[i] = new DataTag();
            dataInfo[i] = new DataFunction.DataInfoFunction("stuff", Null.DIMENSION, new DataSourceIndependentSimple(inDsi.getIndependentData(1).getData(),inDsi.getIndependentDataInfo(1)));
            dataInfo[i].addTag(tags[i]);
            data[i] = new DataFunction(new int[]{shape[1]});
        }
    }
}
