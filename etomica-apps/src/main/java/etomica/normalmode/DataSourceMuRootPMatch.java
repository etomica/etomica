/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.AccumulatorAverageBlockless;
import etomica.data.DataDistributer;
import etomica.data.DataSourceScalar;
import etomica.integrator.mcmove.MCMoveOverlapListener;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure;

/**
 * Datasource that returns an estimate of mu based on thermodynamic
 * self-consistency of pressure measurements and defect free energies.
 * 
 * @author Andrew Schultz
 */
public class DataSourceMuRootPMatch extends DataSourceScalar {

    protected final MCMoveOverlapListener mcMoveOverlapMeter;
    protected final DataDistributer pSplitter;
    protected double bmu;
    protected double bALattice, volume, latticeDensity;
    protected double lastP, lastVacancyConcentration, lastDP, lastLnTot;
    
    public DataSourceMuRootPMatch(MCMoveOverlapListener mcMoveOverlapMeter, double bmu, DataDistributer pSplitter, double bALattice, double latticeDensity, double volume) {
        super("mu", Null.DIMENSION);
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.pSplitter = pSplitter;
        this.bmu = bmu;
        this.latticeDensity = latticeDensity;
        this.bALattice = bALattice;
        this.volume = volume;
    }
    
    public synchronized double getDataAsScalar() {
        double[] ratios = mcMoveOverlapMeter.getRatios();
        if (ratios == null || ratios.length == 0 || pSplitter.getNumDataSinks() < 2) return Double.NaN;
        double myMu = bmu;
        double lastMu = bmu;
        double maxMu = Double.POSITIVE_INFINITY;
        double minMu = Double.NEGATIVE_INFINITY;
        while (true) {

            double p = 1;
            double tot = 0;
            double l = Math.exp(myMu);
            double vAvg = 0;
            double totMinus1 = 0;
            for (int i=ratios.length-1; i>=0; i--) {
                tot += p;
                if (i<ratios.length-1) totMinus1 += p;
                if (Double.isNaN(ratios[i])) {
                    if (i==ratios.length-1) return Double.NaN;
                    break;
                }
                p /= l*ratios[i];
            }
            tot += p;
            if (totMinus1 > 0.5) {
                // probably only happens for large systems where multiple vacanacies are happy
                lastLnTot = Math.log(tot);
            }
            else {
                // ln(1+x) = x - x^2/2 + x^3/3 + ...
                lastLnTot = 0;
                double xn = totMinus1;
                for (int i=1; i<100; i++) {
                    if (Math.abs(xn/i) < lastLnTot*1e-15) break;
                    lastLnTot += xn/i;
                    xn *= -totMinus1;
                }
            }
            p = 1;
            double pressure1 = 0;
            double oneMinusP0 = 0;
            for (int i=0; i<=ratios.length; i++) {
                if (i>0) oneMinusP0 += p/tot;
                if (ratios.length-1-i >= 0) {
                    if (Double.isNaN(ratios[ratios.length-1-i])) {
                        break;
                    }
                    p /= l*ratios[ratios.length-1-i];
                }
            }
            p = 1;
            for (int i=0; i<pSplitter.getNumDataSinks() && i<=ratios.length; i++) {
                double pi = p/tot;
                AccumulatorAverageBlockless acc = (AccumulatorAverageBlockless)pSplitter.getDataSink(i);
                if (acc == null || acc.getSampleCount() == 0) {
                    if (i<2) return Double.NaN;
                    break;
                }
                pressure1 += pi*acc.getData().getValue(acc.AVERAGE.index);
                if (i==0) lastDP = -oneMinusP0*acc.getData().getValue(acc.AVERAGE.index);
                else lastDP += pi*acc.getData().getValue(acc.AVERAGE.index);
                vAvg += pi*i;
                if (ratios.length-1-i >= 0) {
                    if (Double.isNaN(ratios[ratios.length-1-i])) {
                        break;
                    }
                    p /= l*ratios[ratios.length-1-i];
                }
            }
            lastVacancyConcentration = vAvg / (volume*latticeDensity);

            // we could take pLattice from the data collected when nVacancy=0,
            // but for a large system, that might rarely happen
            double pressure2 = (myMu-bALattice)*latticeDensity+Math.log(tot)/volume;
            lastP = pressure2;
            if (Double.isNaN(pressure1) || Double.isNaN(pressure2)) {
                throw new RuntimeException("oops");
            }
            double newMu = myMu - (pressure2-pressure1)/latticeDensity;
//            if (newMu < 15) {
//                System.out.println("oops");
//            }
            if (newMu == myMu || newMu == lastMu || newMu > maxMu || newMu < minMu) {
                return newMu;
            }
            if (((newMu-myMu)*(lastMu-myMu) > 0 && (Math.abs(newMu-myMu) >= Math.abs(lastMu-myMu))) || newMu > maxMu || newMu < minMu) {
                // getting worse
                return myMu;
            }
            if (newMu < myMu) {
                // myMu was too large
                if (maxMu > myMu) maxMu = myMu;
            }
            else {
                // myMu was too small
                if (minMu < myMu) minMu = myMu;
            }
            lastMu = myMu;
            myMu = newMu;
        }
    }
    
    public double getLastPressure() {
        return lastP;
    }
    
    public double getLastDPressure() {
        return lastDP;
    }

    public double getLastVacancyConcentration() {
        return lastVacancyConcentration;
    }
    
    public double getLastLnTot() {
        return lastLnTot;
    }

    public class DataSourceMuRootPressure extends DataSourceScalar {
        public DataSourceMuRootPressure() {
            super("P", Pressure.DIMENSION);
        }
        
        public double getDataAsScalar() {
            DataSourceMuRootPMatch.this.getDataAsScalar();
            return DataSourceMuRootPMatch.this.getLastPressure();
        }
    }

    public class DataSourceMuRootVacancyConcentration extends DataSourceScalar {
        public DataSourceMuRootVacancyConcentration() {
            super("vc", Null.DIMENSION);
        }
        
        public double getDataAsScalar() {
            DataSourceMuRootPMatch.this.getDataAsScalar();
            return DataSourceMuRootPMatch.this.getLastVacancyConcentration();
        }
    }
}
