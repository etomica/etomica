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
 * self-consistency of pressure measurements and defect free energies. This
 * class only considers the free energy of the first defect (it assumes that
 * all additional defects do not interact).  The pressure is also only computed
 * for the perfect and single-vacancy crystal and extrapolated from there.
 * 
 * @author Andrew Schultz
 */
public class DataSourceMuRoot1 extends DataSourceScalar {

    protected final MCMoveOverlapListener mcMoveOverlapMeter;
    protected final DataDistributer pSplitter;
    protected double bmu;
    protected double bALattice, volume, latticeDensity;
    protected double lastP, lastVacancyConcentration;
    protected double minMu;
    
    public DataSourceMuRoot1(MCMoveOverlapListener mcMoveOverlapMeter, double bmu, DataDistributer pSplitter, double bALattice, double latticeDensity, double volume) {
        super("mu", Null.DIMENSION);
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.pSplitter = pSplitter;
        this.bmu = bmu;
        this.latticeDensity = latticeDensity;
        this.bALattice = bALattice;
        this.volume = volume;
        minMu = Double.NaN;
    }
    
    public void setMinMu(double minMu) {
        this.minMu = minMu;
    }
    
    public synchronized double getDataAsScalar() {
        double[] ratios = mcMoveOverlapMeter.getRatios();
        if (ratios == null || ratios.length == 0 || pSplitter.getNumDataSinks() < 2) return Double.NaN;
        double daDef = Math.log(ratios[ratios.length-1]*(volume/latticeDensity));
        if (Double.isNaN(daDef)) return Double.NaN;
        double myMu = bmu;
        double lastMu = bmu;
        double maxMu = Double.POSITIVE_INFINITY;
        double thisMinMu = Double.isNaN(minMu) ? Double.NEGATIVE_INFINITY : minMu;
        int n0 = mcMoveOverlapMeter.getMinNumAtoms();
        int nLatticeAtoms = n0 + ratios.length;
        while (true) {

            double p = 1;
            double tot = 0;
            double vAvg = 0;
            for (int i=0; p/tot > 1e-14; i++) {
                tot += p;
                double x = Math.exp(daDef+myMu)/(nLatticeAtoms-i)*(i+1);
                if (Double.isInfinite(x) || Double.isInfinite(p/x)) throw new RuntimeException("oops "+p+" "+x+" "+daDef+" "+myMu+" "+thisMinMu+" "+i);
                p /= Math.exp(daDef+myMu)/(nLatticeAtoms-i)*(i+1);
                if (Double.isInfinite(p)) throw new RuntimeException("oops");
                if (i>nLatticeAtoms/10) return Double.NaN;
            }
            tot += p;
            p = 1;
            double pressure1 = 0;
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
            double pressure2 = (myMu-bALattice)*latticeDensity+Math.log(tot)/volume;
            lastP = pressure2;
            if (Double.isNaN(pressure2) || Double.isInfinite(pressure2)) {
                throw new RuntimeException("oops "+myMu+" "+bALattice+" "+latticeDensity+" "+tot+" "+volume);
            }
            // calculate a new estimate of mu
            double newMu = myMu - (pressure2-pressure1)/latticeDensity;
            if (Double.isNaN(newMu) || Double.isInfinite(newMu)) {
                throw new RuntimeException("oops "+myMu+" "+pressure1+" "+pressure2+" "+latticeDensity);
            }
//            if (newMu < 15) {
//                System.out.println("oops");
//            }
            if (newMu == myMu || newMu == lastMu || newMu > maxMu || newMu < thisMinMu) {
                return newMu;
            }
            if (((newMu-myMu)*(lastMu-myMu) > 0 && (Math.abs(newMu-myMu) >= Math.abs(lastMu-myMu))) || newMu > maxMu || newMu < thisMinMu) {
                // getting worse
                return myMu;
            }
            if (newMu < myMu) {
                // myMu was too large
                if (maxMu > myMu) maxMu = myMu;
            }
            else {
                // myMu was too small
                if (thisMinMu < myMu) thisMinMu = myMu;
            }
            lastMu = myMu;
            myMu = newMu;
            if (!Double.isNaN(thisMinMu) && myMu < thisMinMu) return Double.NaN;
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
            DataSourceMuRoot1.this.getDataAsScalar();
            return DataSourceMuRoot1.this.getLastPressure();
        }
    }

    public class DataSourceMuRootVacancyConcentration extends DataSourceScalar {
        public DataSourceMuRootVacancyConcentration() {
            super("vc", Null.DIMENSION);
        }
        
        public double getDataAsScalar() {
            DataSourceMuRoot1.this.getDataAsScalar();
            return DataSourceMuRoot1.this.getLastVacancyConcentration();
        }
    }
}
