/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.PotentialMaster;
import etomica.units.dimensions.Null;

/**
 * Meter whose purpose in life is to measure the overlap between systems
 * defined by a combination of two energy meters.  The returned energy is taken
 * to be   (1-f)*meterRef + f*meterTarg
 * The systems under consideration for overlap differ only in the values of f.
 *
 * @author Andrew Schultz
 */
public class MeterOverlapSwitch implements IDataSource {

    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final DataTag tag;
    protected final MeterPotentialEnergy meterRef, meterTarget;
    protected double temperature;
    protected double latticeEnergy;
    protected double frac;
    protected double[] otherFrac;
    protected double[][] alpha;
    protected double[] alphaCenter;
    protected double alphaSpan;
    protected int numAlpha = 1;
    public double targetSum, refSum;
    public int count;
    
    
    public MeterOverlapSwitch(PotentialMaster potentialMasterRef, PotentialMaster potentialMasterTarget) {
        meterRef = new MeterPotentialEnergy(potentialMasterRef);
        meterTarget = new MeterPotentialEnergy(potentialMasterTarget);
        meterTarget.setIncludeLrc(false);
        tag = new DataTag();
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setBox(Box newBox) {
        meterRef.setBox(newBox);
        meterTarget.setBox(newBox);
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public void setSampledSwitchFrac(double newFrac) {
        frac = newFrac;
    }
    
    public void setOtherSwitchFrac(double[] newOtherFrac) {
        otherFrac = newOtherFrac;
    }
    
    protected void initAlpha() {
        if (alphaCenter == null) {
            return;
        }
        alpha = new double[alphaCenter.length][numAlpha];
        for (int i=0; i<alpha.length; i++) {
            if (numAlpha == 1) {
                alpha[i][0] = alphaCenter[i];
            }
            else {
                for (int j=0; j<numAlpha; j++) {
                    alpha[i][j] = alphaCenter[i]*Math.exp(2.0*alphaSpan*(j-(numAlpha-1)/2)/(numAlpha-1));
                }
            }
        }
        data = new DataDoubleArray(numAlpha*alphaCenter.length);
        dataInfo = new DataInfoDoubleArray("overlap", Null.DIMENSION, new int[]{numAlpha*alphaCenter.length});
    }
    
    public double[] getAlpha(int iTemp) {
        return alpha[iTemp];
    }
    
    public void setAlpha(double[] newAlpha) {
        alphaCenter = newAlpha;
        initAlpha();
    }
    
    public void setAlphaSpan(double newAlphaSpan) {
        alphaSpan = newAlphaSpan;
        initAlpha();
    }
    
    public void setNumAlpha(int newNumAlpha) {
        numAlpha = newNumAlpha;
        initAlpha();
    }

    public double getLatticeEnergy() {
        return latticeEnergy;
    }

    public void setLatticeEnergy(double latticeEnergy) {
        this.latticeEnergy = latticeEnergy;
    }

    public IData getData() {
        double uTarget = meterTarget.getDataAsScalar() - latticeEnergy;
        double uRef = meterRef.getDataAsScalar();
        targetSum += uTarget;
        refSum += uRef;
        count++;
        double uSampled = frac*uTarget + (1-frac)*uRef;

        double[] x = data.getData();
        for (int i=0; i<otherFrac.length; i++) {
            double uOther = otherFrac[i]*uTarget + (1-otherFrac[i])*uRef;
            double eRatio = Math.exp(-(uSampled-uOther)/temperature);
            for (int j=0; j<numAlpha; j++) {
                if (frac>otherFrac[i]) {
                    // eOther / (alpha * eOther + eSampled)
                    // 1.0 / (alpha + eSampled/eOther)
                    x[i*numAlpha+j] = 1.0/(alpha[i][j]+eRatio);
                }
                else {
                    x[i*numAlpha+j] = 1.0/(1+alpha[i][j]*eRatio);
                }
            }
        }
        return data;
    }
}
