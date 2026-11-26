/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataDistributer;
import etomica.data.DataSourceScalar;
import etomica.integrator.mcmove.MCMoveOverlapListener;
import etomica.units.dimensions.Pressure;

/**
 * Averages the measured pressures for various N.  The average is weighted by
 * the probability of visiting that N based on the imposed chemical potential
 * and the measured free energy differences.
 * 
 * @author Andrew Schultz
 */
public class DataSourceAvgPressure extends DataSourceScalar {
    protected final DataDistributer pSplitter;
    protected final MCMoveOverlapListener mcMoveOverlapMeter;
    protected double bmu;

    public DataSourceAvgPressure(DataDistributer pSplitter, MCMoveOverlapListener mcMoveOverlapMeter, double bmu) {
        super("pressure", Pressure.DIMENSION);
        this.pSplitter = pSplitter;
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.bmu = bmu;
    }

    public void setMu(double newBMu) {
        bmu = newBMu;
    }

    public double getDataAsScalar() {
        if (pSplitter.getNumDataSinks() == 0 || pSplitter.getDataSink(0) == null || ((AccumulatorAverageCollapsing)pSplitter.getDataSink(0)).getSampleCount() == 0) return Double.NaN;
        double[] ratios = mcMoveOverlapMeter.getRatios();
        if (ratios == null) return Double.NaN;
        double p = 1;
        double tot = 0;
        for (int i=ratios.length-1; i>=0; i--) {
            tot += p;
            if (Double.isNaN(ratios[i])) {
                break;
            }
            p *= Math.exp(-bmu)/ratios[i];
        }
        tot += p;
        double p2 = 1;
        double pressure = 0;
        for (int i=0; i<pSplitter.getNumDataSinks() && i<=ratios.length; i++) {
            double pi = p2/tot;
            AccumulatorAverageCollapsing acc = (AccumulatorAverageCollapsing)pSplitter.getDataSink(i);
            if (acc == null || acc.getSampleCount() == 0) {
                break;
            }
            pressure += pi*acc.getData().getValue(acc.AVERAGE.index);
            if (ratios.length-1-i >= 0) {
                p2 *= Math.exp(-bmu)/ratios[ratios.length-1-i];
                if (Double.isNaN(ratios[ratios.length-1-i])) {
                    break;
                }
            }
        }
        return pressure;
    }
}
