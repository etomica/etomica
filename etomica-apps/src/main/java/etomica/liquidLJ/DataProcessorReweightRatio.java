/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

/**
 * Takes incoming data sums including a weight sum for each cutoff,
 * computes averages as sum/(weight sum)
 *
 * @author andrew
 */
public class DataProcessorReweightRatio extends DataProcessor {

    protected DataDoubleArray data;
    protected final int nCutoffs;
    protected final int ref;
    
    public DataProcessorReweightRatio(int nCutoffs) {
        this(nCutoffs, -1);
    }
    
    public DataProcessorReweightRatio(int nCutoffs, int ref) {
        this.nCutoffs = nCutoffs;
        this.ref = ref;
    }

    protected IData processData(IData inputData) {
        double[] x = data.getData();
        int j = 0;
        int nData = inputData.getLength()/nCutoffs-1;
        for (int i=0; i<nCutoffs; i++) {
            double w = inputData.getValue(j+i+nData);
            for (int k=0; k<nData; k++) {
                x[j+k] = inputData.getValue(j+i+k)/w;
            }
            j += nData;
        }
        if (ref>=0) {
            j=0;
            for (int i=0; i<nCutoffs; i++) {
                if (i!=ref) {
                    for (int k=0; k<nData; k++) {
                        x[j+k] -= x[ref*nData+k];
                    }
                }
                j += nData;
           }
        }
        if (data.isNaN()) {
            throw new RuntimeException("oops");
        }
        return data;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        int nData = inputDataInfo.getLength()/nCutoffs-1;
        dataInfo = new DataInfoDoubleArray("whatever", Null.DIMENSION, new int[]{nData*nCutoffs});
        data = new DataDoubleArray(dataInfo.getLength());
        return dataInfo;
    }
}
