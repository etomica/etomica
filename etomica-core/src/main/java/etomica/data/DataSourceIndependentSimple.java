/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;

/**
 * A simple implementation of DataSourceIndependent objects can use if they 
 * cannot act as the DataSourceIndependent themselves (if perhaps they make
 * multiple DataFunctions)
 */
public class DataSourceIndependentSimple implements DataSourceIndependent {

    protected final DataTag tag;
    private DataDoubleArray[] xData;
    private DataInfoDoubleArray[] xDataInfo;

    public DataSourceIndependentSimple(double[] rawData, DataInfoDoubleArray xDataInfo) {
        this(new double[][]{rawData}, new DataInfoDoubleArray[]{xDataInfo});
    }

    public DataSourceIndependentSimple(double[][] rawData, DataInfoDoubleArray[] xDataInfo) {
        tag = new DataTag();
        xData = new DataDoubleArray[rawData.length];
        this.xDataInfo = xDataInfo;
        for (int i=0; i<rawData.length; i++) {
            xData[i] = new DataDoubleArray(new int[]{rawData[i].length}, rawData[i]);
            xDataInfo[i].addTag(tag);
        }
    }
    
    public void update(double[][] rawData, DataInfoDoubleArray[] newXDataInfo) {
        System.arraycopy(xDataInfo, 0, xDataInfo, 0, newXDataInfo.length);
        for (int i=0; i<rawData.length; i++) {
            xData[i] = new DataDoubleArray(new int[]{rawData[i].length}, rawData[i]);
            xDataInfo[i].addTag(tag);
        }
    }

    public void update(double[] rawData, DataInfoDoubleArray newXDataInfo) {
        this.xDataInfo[0] = newXDataInfo;
        xData[0] = new DataDoubleArray(new int[]{rawData.length}, rawData);
        this.xDataInfo[0].addTag(tag);
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo[i];
    }
    
    public DataDoubleArray getIndependentData(int i) {
        return xData[i];
    }

    public int getIndependentArrayDimension() {
        return xData.length;
    }

    public DataTag getIndependentTag() {
        return tag;
    }
}