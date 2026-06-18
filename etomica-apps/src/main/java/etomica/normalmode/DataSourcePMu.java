/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.integrator.mcmove.MCMoveOverlapListener;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure;

/**
 * DataSource that returns P2-P1 or mu-mu* as a function of mu.  For either
 * case, 0 means thermodynamic self-consistency.
 * 
 * P1 = <P>
 * P2 = (mu-beta*ALattice)*latticeDensity+Math.log(sum)/volume;
 *   sum = sum[Pi/Pn], Pi=probability of visiting N=i, n=# lattice site
 * mu* = (P2-P1)/latticeDensity
 * 
 * @author Andrew Schultz
 */
public class DataSourcePMu implements IDataSource, DataSourceIndependent {

    protected final MCMoveOverlapListener mcMoveOverlapMeter;
    protected final DataDistributer pSplitter;
    protected DataFunction data;
    protected DataInfoFunction dataInfo;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    protected final DataTag tag, xTag;
    protected double bmu, deltaBMu;
    protected int nMu;
    protected double bALattice, volume, latticeDensity;
    protected boolean returnMu;

    public DataSourcePMu(MCMoveOverlapListener mcMoveOverlapMeter, double deltaBMu, int nMu, double bmu, DataDistributer pSplitter, double bALattice, double latticeDensity, double volume, boolean returnMu) {
        tag = new DataTag();
        xTag = new DataTag();
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.pSplitter = pSplitter;
        xDataInfo = new DataInfoDoubleArray("betaMu", Null.DIMENSION, new int[]{nMu});
        xData = new DataDoubleArray(nMu);
        this.deltaBMu = deltaBMu;
        this.nMu = nMu;
        setMu(bmu);
        data = new DataFunction(new int[]{nMu});
        dataInfo = new DataInfoFunction("deltaP", Pressure.DIMENSION, this);
        this.latticeDensity = latticeDensity;
        this.bALattice = bALattice;
        this.volume = volume;
        this.returnMu = returnMu;
    }
    
    public void setMu(double newMu) {
        double[] x = xData.getData();
        bmu = newMu;
        for (int i=-nMu/2; i<=nMu/2; i++) {
            x[i+nMu/2] = bmu+i*deltaBMu;
        }
        xDataInfo = new DataInfoDoubleArray("betaMu", Null.DIMENSION, new int[]{nMu});
    }
    
    public IData getData() {
        double[] ratios = mcMoveOverlapMeter.getRatios();
        if (ratios == null) return data;
        double[] y = data.getData();
        for (int j=0; j<nMu; j++) {
            double imu = xData.getData()[j];

            double p = 1;
            double tot = 0;
            double l = Math.exp(imu);
            for (int i=ratios.length-1; i>=0; i--) {
                tot += p;
                if (Double.isNaN(ratios[i])) {
                    break;
                }
                p *= 1/(l*ratios[i]);
            }
            tot += p;
            double p2 = 1;
            double pressure1 = 0;
            for (int i=0; i<pSplitter.getNumDataSinks() && i<=ratios.length; i++) {
                double pi = p2/tot;
                AccumulatorAverageBlockless acc = (AccumulatorAverageBlockless)pSplitter.getDataSink(i);
                if (acc == null || acc.getSampleCount() == 0) {
                    if (i==0) return data;
                    break;
                }
                pressure1 += pi*acc.getData().getValue(acc.AVERAGE.index);
                if (ratios.length-1-i >= 0) {
                    p2 *= 1/(l*ratios[ratios.length-1-i]);
                    if (Double.isNaN(ratios[ratios.length-1-i])) {
                        break;
                    }
                }
            }

            double q = 1;
            double sum = 1;
            for (int i=ratios.length-1; i>=0; i--) {
                if (Double.isNaN(ratios[i])) break;
                q *= 1/(ratios[i]*l);
                sum += q;
            }
            // we could take pLattice from the data collected when nVacancy=0,
            // but for a large system, that might rarely happen
            double pressure2 = (imu-bALattice)*latticeDensity+Math.log(sum)/volume;
            if (returnMu) {
                y[j] = imu - (pressure2-pressure1)/latticeDensity;
            }
            else {
                y[j] = pressure2-pressure1;
            }
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
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
