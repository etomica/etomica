/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.box.Box;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

/**
 * Takes incoming data, computes weight for configuration based on various
 * cutoffs, multiplies each incoming data by the weight and also adds
 * weights as extra data.
 * 
 * Taking n incoming data and m cutoffs, (n+1)*m data will be sent down the
 * pipe.
 * 
 * @author andrew
 */
public class DataProcessorReweight extends DataProcessor {
    private final double temperature;
    private final ValueCache energyFastCache;
    private final double[] uFac;
    protected DataDoubleArray data;
    protected final Box box;
    protected final int nCutoffs;

    public DataProcessorReweight(double temperature, ValueCache energyFastCache,
                                 double[] uFac, Box box, int nCutoffs) {
        this.temperature = temperature;
        this.energyFastCache = energyFastCache;
        this.uFac = uFac;
        this.box = box;
        this.nCutoffs = nCutoffs;
    }

    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = new DataInfoDoubleArray("whatever", Null.DIMENSION, new int[]{(inputDataInfo.getLength()+nCutoffs)});
        data = new DataDoubleArray(dataInfo.getLength());
        return dataInfo;
    }

    protected IData processData(IData inputData) {
        double uFast = energyFastCache.getValue();
        double[] x = data.getData();
        int j = 0;
        int nData = inputData.getLength()/nCutoffs;
        int n = box.getMoleculeList().getMoleculeCount();
        for (int i=0; i<nCutoffs; i++) {
            double dx = n*inputData.getValue(j) - (uFast+uFac[i]);
//            System.out.println(i+" "+n*inputData.getValue(j)+" "+uFast+" "+uFac[i]);
            double w = Math.exp(-dx/temperature);
            for (int k=0; k<nData; k++) {
                x[j+i+k] = inputData.getValue(j+k)*w;
            }
            x[j+i+nData] = w;
            j += inputData.getLength()/nCutoffs;
        }
//        System.out.println("rw "+Arrays.toString(x));
        if (data.isNaN()) {
            throw new RuntimeException("oops");
        }
        return data;
    }
}
