/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic.mbar;

import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

public class MeterMBAR implements IEtomicaDataSource, IntegratorListener {

    protected double[] alphaSpan, alphaCenter;
    protected int nAlpha;
    protected double[][] alpha;
    protected final DataTag tag;
    protected DataDoubleArray data;
    protected double[] sumCheck;
    protected DataInfoDoubleArray dataInfo;
    protected final DataSourceScalar[] meterPE;
    protected double temperature;
    protected final double[] ei;
    protected int index;
    protected long callCount;
    
    public MeterMBAR(DataSourceScalar[] meterPE, double temperature) {
        this.meterPE = meterPE;
        this.temperature = temperature;
        tag = new DataTag();
        ei = new double[2];
    }

    public int getNumAlpha() {
        return nAlpha;
    }
    
    public double[] getAlphaCenter() {
        return alphaCenter;
    }
    
    public double[] getAlphaSpan() {
        return alphaSpan;
    }

    public double[] getAlpha(int iAlpha) {
        return alpha[iAlpha];
    }
    
    public long getCallCount() {
        return callCount;
    }

    public void reset() {
        data.E(0);
        callCount = 0;
    }
    
    public void setNumAlpha(int newNumAlpha) {
        nAlpha = newNumAlpha;
        data = new DataDoubleArray(new int[]{nAlpha,meterPE.length});
        dataInfo = new DataInfoDoubleArray("chi", Null.DIMENSION, new int[]{nAlpha,meterPE.length});
        alpha = new double[nAlpha][meterPE.length-1];
        if (alphaCenter != null) {
            setAlpha(alphaCenter, alphaSpan);
        }
        sumCheck = new double[nAlpha];
    }
    
    /**
     * sets the range of parameter values used for Bennets method.
     * Default is a span of 5 centered about 1 (exp(-5) to (exp(5)).
     * @param aCenter geometric mean of all values
     * @param aSpan natural log of ratio of max value to aCenter
     */
    public void setAlpha(double[] aCenter, double[] aSpan) {
        for (int i=0; i<aCenter.length; i++) {
            if (aSpan[i] < 0.0 || (aSpan[i] == 0 && nAlpha > 1) || aCenter[i] <= 0.0 ) throw new IllegalArgumentException("span and center must be positive");
        }
        alphaSpan = aSpan;
        alphaCenter = aCenter;
        if (nAlpha == 0) {
            return;
        }
        if (nAlpha==1) {
            alpha[0] = aCenter;
        }
        else {
            for (int i=0; i<aCenter.length; i++) {
                for (int j=0; j<nAlpha; j++) {
                    alpha[j][i] = Math.exp(2.0*aSpan[i]*(j-(nAlpha-1)/2)/(nAlpha-1))*aCenter[i];
                }
            }
        }
        reset();
    }
    
    public IData getData() {
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public void setBoxIndex(int newIndex) {
        index = newIndex;
    }
    
    public int getBoxIndex() {
        return index;
    }

    public void integratorStepFinished(IntegratorEvent e) {
        callCount++;
        double[] sum = data.getData();
        int n = meterPE.length;
        for (int j=0; j<n; j++) {
            ei[j] = Math.exp(-meterPE[j].getDataAsScalar()/temperature);
        }
        for (int i=0; i<nAlpha; i++) {
            sumCheck[i] += ei[1-index] / (ei[0] + alpha[i][0]*ei[1]);
            
            double denom = ei[0];
            for (int j=1; j<n; j++) {
                denom += alpha[i][j-1]*ei[j];
            }
            double total = 0;
            for (int j=0; j<n; j++) {
                if (j != index) {
                    sum[2*i+j] += ei[j]/denom;
                    total += ei[j];
                }
            }
            if (index == 0) {
                sum[n*i+index] -= total/denom;
            }
            else {
                sum[n*i+index] -= total/(denom*alpha[i][index-1]);
            }
            
//            System.out.println("check "+index+" "+i+" "+sum[2*i+(1-index)]+" "+sumCheck[i]);
        }
    }

    public void integratorInitialized(IntegratorEvent e) {}
    public void integratorStepStarted(IntegratorEvent e) {}
}
