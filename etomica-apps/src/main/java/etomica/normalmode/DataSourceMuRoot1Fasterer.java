/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.AccumulatorAverageBlockless;
import etomica.data.DataDistributer;
import etomica.data.DataSourceScalar;
import etomica.integrator.mcmove.MCMoveOverlapListenerFasterer;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure;

/**
 * Datasource that returns an estimate of mu based on thermodynamic
 * self-consistency of pressure measurements and defect free energies. This
 * class only considers the free energy of the first defect (it assumes that
 * all additional defects do not interact).  The pressure is also only computed
 * for the perfect and single-vacancy crystal and extrapolated from there.
 * 
 * @author Andrew Schultz
 */
public class DataSourceMuRoot1Fasterer extends DataSourceScalar {

    protected final MCMoveOverlapListenerFasterer mcMoveOverlapMeter;
    protected final DataDistributer pSplitter;
    protected double bmu;
    protected double bALattice, volume, latticeDensity;
    protected double lastP, lastVacancyConcentration;
    protected double lastDP, lastLnTot;

    public DataSourceMuRoot1Fasterer(MCMoveOverlapListenerFasterer mcMoveOverlapMeter, double bmu, DataDistributer pSplitter, double bALattice, double latticeDensity, double volume) {
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
        double daDef = Math.log(ratios[ratios.length-1]*(volume/latticeDensity));
        if (Double.isNaN(daDef)) return Double.NaN;
        double myMu = bmu;
        double lastMu = bmu;
        double maxMu = Double.POSITIVE_INFINITY;
        double minMu = Double.NEGATIVE_INFINITY;
        double thisMinMu = Double.isNaN(minMu) ? Double.NEGATIVE_INFINITY : minMu;
        int n0 = mcMoveOverlapMeter.getMinNumAtoms();
        int nLatticeAtoms = n0 + ratios.length;
        int iters = 0;
        while (true) {

            double p = 1;
            double totMinus1 = 0;
            double vAvg = 0;
            for (int i = 0; p > 1e-14; i++) {
                double x = Math.exp(daDef+myMu)/(nLatticeAtoms-i)*(i+1);
                if (Double.isInfinite(x) || Double.isInfinite(p/x)) throw new RuntimeException("oops "+p+" "+x+" "+daDef+" "+myMu+" "+thisMinMu+" "+i);
                p /= Math.exp(daDef+myMu)/(nLatticeAtoms-i)*(i+1);
                if (Double.isInfinite(p)) throw new RuntimeException("oops");
                if (i>nLatticeAtoms/10) return Double.NaN;
                totMinus1 += p;
            }

            double tot = totMinus1 + 1;
            if (totMinus1 > 0.5) {
                // probably only happens for large systems where multiple vacanacies are happy
                lastLnTot = Math.log(tot);
            } else {
                // ln(1+x) = x - x^2/2 + x^3/3 + ...
                lastLnTot = 0;
                double xn = totMinus1;
                for (int i = 1; i < 100; i++) {
                    if (Math.abs(xn / i) < lastLnTot * 1e-15) break;
                    lastLnTot += xn / i;
                    xn *= -totMinus1;
                }
            }

            p = 1;
            double pressure1 = 0;
            double oneMinusP0 = 0;
            for (int i = 0; i <= ratios.length; i++) {
                if (i > 0) oneMinusP0 += p / tot;
                p /= Math.exp(daDef + myMu) / (nLatticeAtoms - i) * (i + 1);
            }
            p = 1;
            AccumulatorAverageBlockless acc0 = (AccumulatorAverageBlockless)pSplitter.getDataSink(0);
            AccumulatorAverageBlockless acc1 = (AccumulatorAverageBlockless)pSplitter.getDataSink(1);
            if (acc0 == null || acc1 == null || acc0.getSampleCount() == 0 || acc1.getSampleCount() == 0) return Double.NaN;
            double p0 = acc0.getData().getValue(acc0.AVERAGE.index);
            double p1 = acc1.getData().getValue(acc1.AVERAGE.index);
            double dp = p1-p0;
            if (dp > 0) dp = 0;
            for (int i=0; p/tot > 1e-14; i++) {
                double pi = p/tot;
                double iPressure = p0 + dp*i;
                if (i == 0) {
                    lastDP = -oneMinusP0 * iPressure;
                } else {
                    lastDP += pi * iPressure;
                }
                pressure1 += pi*iPressure;
                if (Double.isNaN(pressure1)) {
                    throw new RuntimeException("oops "+pi+" "+iPressure);
                }
                vAvg += pi*i;
                p /= Math.exp(daDef+myMu)/(nLatticeAtoms-i)*(i+1);
            }
            lastVacancyConcentration = vAvg / (volume*latticeDensity);

            // we could take pLattice from the data collected when nVacancy=0,
            // but for a large system, that might rarely happen
            double dmu = lastDP / latticeDensity - lastLnTot / volume;
            lastP = pressure1;
            // calculate a new estimate of mu
            double newMu = bmu + dmu;
            if (Double.isNaN(newMu) || Double.isInfinite(newMu)) {
                throw new RuntimeException("oops " + myMu + " " + lastDP + " " + dmu + " " + latticeDensity);
            }
            if (newMu == myMu || newMu == lastMu) {
                return newMu;
            }
            if (newMu >= maxMu) {
                newMu = (myMu + maxMu) / 2;
                if (newMu == maxMu || newMu == myMu) return newMu;
            }
            if (newMu <= minMu) {
                newMu = (minMu + myMu) / 2;
                if (newMu == minMu || newMu == myMu) return newMu;
            }
            if (newMu > myMu) {
                minMu = Math.max(myMu, minMu);
            }
            else if (newMu < myMu){
                maxMu = Math.min(myMu, maxMu);
            }
            if ((newMu - myMu) * (lastMu - myMu) > 0 && iters > 10) {
                // we're bouncing.  if bouncing hard, consider going to the middle
                if ((newMu - myMu) / (lastMu - myMu) > 0.9) {
                    newMu = (newMu + myMu) / 2;
                }
            }
            lastMu = myMu;
            myMu = newMu;
            iters++;
            if (iters>1000) throw new RuntimeException("oops");
        }
    }
    
    public double getLastPressure() {
        return lastP;
    }
    
    public double getLastVacancyConcentration() {
        return lastVacancyConcentration;
    }

    public class DataSourceMuRootPressure extends DataSourceScalar {
        public DataSourceMuRootPressure() {
            super("P", Pressure.DIMENSION);
        }

        public double getDataAsScalar() {
            DataSourceMuRoot1Fasterer.this.getDataAsScalar();
            return DataSourceMuRoot1Fasterer.this.getLastPressure();
        }
    }

    public class DataSourceMuRootVacancyConcentration extends DataSourceScalar {
        public DataSourceMuRootVacancyConcentration() {
            super("vc", Null.DIMENSION);
        }
        
        public double getDataAsScalar() {
            DataSourceMuRoot1Fasterer.this.getDataAsScalar();
            return DataSourceMuRoot1Fasterer.this.getLastVacancyConcentration();
        }
    }
}
