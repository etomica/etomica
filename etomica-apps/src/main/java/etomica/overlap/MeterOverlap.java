/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.overlap;

import etomica.data.DataTag;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;


public class MeterOverlap implements IDataSource, AlphaSource {
    protected DataTag tag;
    protected IEtomicaDataInfo dataInfo;
    protected IDataSource dataSourceRef, dataSourceTarget;
    protected double temperature;
    protected DataDoubleArray data;
    protected int numAlpha;
    protected double alphaCenter, alphaSpan;
    protected double[] alpha;
    protected boolean isReference;
    
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
    public MeterOverlap(IDataSource dataSourceRef,
                        IDataSource dataSourceTarget, double temperature, boolean isReference) {
        tag = new DataTag();
        this.dataSourceRef = dataSourceRef;
        this.dataSourceTarget = dataSourceTarget;
        this.temperature = temperature;
        this.isReference = isReference;
        setNumAlpha(5);
        setAlphaRange(1, 5);
    }

    /**
     * sets the range of parameter values used for Bennets method.
     * @param aCenter geometric mean of all values
     * @param aSpan natural log of ratio of max value to aCenter
     */
    public void setAlphaRange(double aCenter, double aSpan) {
        if (aSpan < 0.0 || aSpan == 0 || aCenter <= 0.0 ) throw new IllegalArgumentException("span and center must be positive");
        alphaSpan = aSpan;
        alphaCenter = aCenter;
        if (numAlpha==1) {
            alpha[0] = aCenter;
        }
        else {
            for (int i=0; i<numAlpha; i++) {
                alpha[i] = Math.exp(2.0*aSpan*(i-(numAlpha-1)/2)/(numAlpha-1))*aCenter;
            }
        }

        // make new dataInfo so downstream knows we've changed
        dataInfo = new DataInfoDoubleArray("overlap", Null.DIMENSION, new int[]{numAlpha});
        dataInfo.addTag(tag);
    }
    
    public void setNumAlpha(int newNumAlpha) {
        if (newNumAlpha < 1) {
            throw new RuntimeException("must be positive");
        }
        numAlpha = newNumAlpha;
        alpha = new double[numAlpha];
        if (alphaCenter > 0) {
            setAlphaRange(alphaCenter, alphaSpan);
        }
        else {
            // make new dataInfo so downstream knows we've changed
            dataInfo = new DataInfoDoubleArray("overlap", Null.DIMENSION, new int[]{numAlpha});
            dataInfo.addTag(tag);
        }
        data = new DataDoubleArray(numAlpha);
    }
    
    public double getAlphaCenter() {
        return alphaCenter;
    }
    
    public double getAlphaSpan() {
        return alphaSpan;
    }

    public int getNumAlpha() {
        return numAlpha;
    }

    public double getAlpha(int iAlpha) {
        return alpha[iAlpha];
    }

    public DataDoubleArray getData(){
        // the real work.  calculate overlap function for each alpha
        
        double eRef = Math.exp(-dataSourceRef.getData().getValue(0)/temperature);
        double eTarget = Math.exp(-dataSourceTarget.getData().getValue(0)/temperature);
        double[] x = data.getData();
        for (int i=0; i<numAlpha; i++) {
            x[i] = isReference ? eTarget : eRef;
            x[i] /= (alpha[i] * eRef + eTarget);
            if (x[i] == 0) {
                // we can't take ln(0), so take x=MIN_VALUE insetad.
                // this won't cause real problems because it will go away when we sum it up,
                //   foo + MIN_VALUE = foo
                x[i] = Double.MIN_VALUE;
                // if Double.MIN_VALUE is actually a value that can make a significant contribution
                // then your calculation is almost certainly broken anyway because you need a data type that can handle
                // smaller numbers.
            }
        }
        return data;
    }
    
    public IEtomicaDataInfo getDataInfo(){
        return dataInfo;
    }
    public DataTag getTag() {
        return tag;
    }
}
