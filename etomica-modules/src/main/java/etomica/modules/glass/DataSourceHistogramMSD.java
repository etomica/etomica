/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.data.*;
import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.integrator.IntegratorMD;
import etomica.units.dimensions.Null;

public class DataSourceHistogramMSD implements IDataSource, DataSourceIndependent, DataSourceMSD.MSDSink {

    protected DataDoubleArray xData;
    protected DataDoubleArray.DataInfoDoubleArray xDataInfo;
    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected final DataTag xTag, tag;
    protected final IntegratorMD integrator;
    protected long step0;
    protected boolean enabled;
    protected int interval;
    protected Histogram[] histograms;

    public DataSourceHistogramMSD(IntegratorMD integrator) {
        this.integrator = integrator;
        histograms = new Histogram[60];
        for (int i = 0; i < histograms.length; i++) {
            histograms[i] = new HistogramCollapsing(100);
        }
        tag = new DataTag();
        xTag = new DataTag();
        resetStep0();
    }

    public void resetStep0() {
        step0 = integrator.getStepCount();
        setInterval(0);
    }

    public void setInterval(int i) {
        interval = i;
        double[] y = histograms[i].getHistogram();
        data = new DataFunction(new int[]{y.length}, y);
        double[] x = histograms[i].xValues();
        xData = new DataDoubleArray(new int[]{x.length}, x);
        xDataInfo = new DataDoubleArray.DataInfoDoubleArray("MSD", Null.DIMENSION, new int[]{x.length});
        xDataInfo.addTag(xTag);
        dataInfo = new DataFunction.DataInfoFunction("MSD histogram", Null.DIMENSION, this);
        dataInfo.addTag(tag);
    }

    public void setEnabled(boolean isEnabled) {
        enabled = isEnabled;
        resetStep0();
    }

    public boolean getEnabled() {
        return enabled;
    }

    @Override
    public IData getData() {
        histograms[interval].getHistogram();
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return xTag;
    }

    public void putMSD(int log2interval, long step, double msd) {
        histograms[log2interval].addValue(msd);
    }
}
