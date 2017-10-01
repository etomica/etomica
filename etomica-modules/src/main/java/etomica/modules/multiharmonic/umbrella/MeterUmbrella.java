/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.umbrella;

import etomica.data.DataTag;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;


public class MeterUmbrella implements IDataSource {
    protected DataTag tag;
    protected IEtomicaDataInfo dataInfo;
    protected IDataSource dataSourceA, dataSourceB;
    protected double temperature;
    protected DataDoubleArray dda;
    protected double alpha, alphaFac;
    
    /**
     * Meter to measure the "overlap" function for free energy calculations,
     * used as input to an AccumulatorVirialOverlapSingleAverage.
     * 
     * @param dataSourceSame - data source that returns the energy of the
     *    system being sampled
     * @param dataSourceDifferent - data source that returns the energy of the
     *    system not being sampled
     * @param temperature - the temperature
     */
    public MeterUmbrella(IDataSource dataSourceSame,
                         IDataSource dataSourceDifferent, double temperature) {
        dataInfo = new DataInfoDoubleArray("umbrella", Null.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        
        this.dataSourceA = dataSourceSame;
        this.dataSourceB = dataSourceDifferent;
        this.temperature = temperature;
        
        dda = new DataDoubleArray(2);
        alphaFac = 1;
    }

    public double getAlphaFac() {
        return alphaFac;
    }
    
    public void setAlphaFac(double newAlphaFac) {
        alphaFac = newAlphaFac;
    }
    
    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double newAlpha) {
        alpha = newAlpha;
    }
    
    public DataDoubleArray getData(){
        double[] x = dda.getData();

        double eA = Math.exp(-dataSourceA.getData().getValue(0)/temperature);
        double eB = Math.exp(-dataSourceB.getData().getValue(0)/temperature);
        double pi = eB + alphaFac*alpha*eA;
        x[0] = eA / pi;
        x[1] = eB / pi;
        return dda;
    }
    
    public IEtomicaDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
}
