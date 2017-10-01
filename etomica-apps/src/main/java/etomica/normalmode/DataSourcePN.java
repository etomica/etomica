/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure;
import etomica.units.dimensions.Quantity;

/**
 * DataSource that returns measured pressure as a function of N
 * 
 * @author Andrew Schultz
 */
public class DataSourcePN implements IDataSource, DataSourceIndependent {

    protected final DataDistributer pSplitter;
    protected DataFunction data;
    protected DataInfoFunction dataInfo;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    protected final DataTag tag, xTag;
    protected final int nominalNumAtoms;
    
    public DataSourcePN(DataDistributer pSplitter, int nominalNumAtoms) {
        this.nominalNumAtoms = nominalNumAtoms;
        tag = new DataTag();
        xTag = new DataTag();
        this.pSplitter = pSplitter;
        xDataInfo = new DataInfoDoubleArray("N", Quantity.DIMENSION, new int[]{0});
        xData = new DataDoubleArray(0);
        data = new DataFunction(new int[]{0});
        dataInfo = new DataInfoFunction("pressure", Pressure.DIMENSION, this);
    }
    
    public IData getData() {
        if (pSplitter.getNumDataSinks() == 0) return data;
        if (pSplitter.getNumDataSinks() > dataInfo.getLength()) {
            getDataInfo();
        }
        double[] y = data.getData();
        for (int j=0; j<pSplitter.getNumDataSinks(); j++) {
            AccumulatorAverageBlockless acc = (AccumulatorAverageBlockless)pSplitter.getDataSink(j);
            if (acc == null || acc.getSampleCount() == 0) {
                y[j] = Double.NaN;
                continue;
            }
            y[j] = acc.getData().getValue(acc.AVERAGE.index);
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        if (dataInfo != null && pSplitter.getNumDataSinks() > dataInfo.getLength()) {
            xDataInfo = new DataInfoDoubleArray("N", Quantity.DIMENSION, new int[]{pSplitter.getNumDataSinks()});
            xDataInfo.addTag(tag);
            xData = new DataDoubleArray(xDataInfo.getLength());
            double[] x = xData.getData();
            for (int j=0; j<xDataInfo.getLength(); j++) {
                x[j] = nominalNumAtoms-j;
            }
            
            dataInfo = new DataInfoFunction("pressure", Null.DIMENSION, this);
            dataInfo.addTag(tag);
            data = new DataFunction(new int[]{xDataInfo.getLength()});
        }
        return dataInfo;
    }

    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return xTag;
    }

}
