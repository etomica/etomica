/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Quantity;

/**
 * Calculates dy/dlnx vs. x for the incoming DataFunction.  The derivatives are
 * calculated as y'=(y4-y3)/(ln(x4)-ln(x3)), with x'=sqrt(x3*x4)
 *
 * @author Andrew Schultz
 */
public class DataProcessorDyDLnx extends DataProcessor implements DataSourceIndependent {

    public DataProcessorDyDLnx() {
        nTag = new DataTag();
    }

    public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
        return null;
    }

    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        nDataSource = ((DataInfoFunction)inputDataInfo).getXDataSource();
        int myLength = inputDataInfo.getLength()-1;
        if (myLength < 0) {
            myLength = 0;
        }
        if (outNData != null && myLength == outNData.getLength()) {
            return dataInfo;
        }
        outNData = new DataDoubleArray(myLength);
        outNDataInfo = new DataInfoDoubleArray("block size", Quantity.DIMENSION, new int[]{myLength});
        outNDataInfo.addTag(nTag);
        IData nData = nDataSource.getIndependentData(0);
        double[] x = outNData.getData();
        for (int j=0; j<x.length; j++) {
            x[j] = Math.sqrt(nData.getValue(j)*nData.getValue(j+1));
        }
        dataInfo = new DataInfoFunction("dydlnx", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        data = new DataFunction(new int[]{myLength});
        return dataInfo;
    }

    protected IData processData(IData inputData) {
        if (inputData.getLength() < 1) {
            return data;
        }

        IData nData = nDataSource.getIndependentData(0);
        double[] d = data.getData();
        for (int i=0; i<d.length; i++) {
            double lnN1 = Math.log(nData.getValue(i));
            double lnN2 = Math.log(nData.getValue(i+1));
            d[i] = (inputData.getValue(i+1) - inputData.getValue(i))/(lnN2-lnN1);
        }
        return data;
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataDoubleArray getIndependentData(int i) {
        return outNData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return outNDataInfo;
    }

    public DataTag getIndependentTag() {
        return nTag;
    }

    private static final long serialVersionUID = 1L;
    protected DataFunction data;
    protected DataInfoFunction dataInfo;
    protected DataDoubleArray outNData;
    protected DataInfoDoubleArray outNDataInfo;
    protected DataSourceIndependent nDataSource;
    protected final DataTag nTag;
}
