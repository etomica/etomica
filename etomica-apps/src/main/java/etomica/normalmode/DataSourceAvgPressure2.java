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
 * Computes pressure from the grand canonical route based on measured free
 * energy differences and simulation parameters.
 * 
 * @author Andrew Schultz
 */
public class DataSourceAvgPressure2 extends DataSourceScalar {
    protected final MCMoveOverlapListener mcMoveOverlapMeter;
    protected final DataDistributer pSplitter;
    protected double bmu;
    protected double latticeDensity;
    protected double volume;
    protected DataSourceMuRoot dsmr;

    public DataSourceAvgPressure2(DataDistributer pSplitter, MCMoveOverlapListener mcMoveOverlapMeter, double latticeDensity, double volume) {
        super("pressure", Pressure.DIMENSION);
        this.pSplitter = pSplitter;
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.latticeDensity = latticeDensity;
        this.volume = volume;
    }

    public void setDataSourceMuRoot(DataSourceMuRoot dsmr) {
    	this.dsmr = dsmr;
    }

    public double getDataAsScalar() {
        double[] ratios = mcMoveOverlapMeter.getRatios();
        if (ratios == null) return Double.NaN;
        double q = 1;
        double sum = 1;
        double l = Math.exp(bmu);
        for (int i=ratios.length-1; i>=0; i--) {
            q *= ratios[i]/l;
            sum += q;
        }
        double lastDMu = dsmr.getLastDMu();
        AccumulatorAverageCollapsing acc = (AccumulatorAverageCollapsing)pSplitter.getDataSink(0);
        if (acc == null || acc.getSampleCount() == 0) {
            return Double.NaN;
        }
        double pLat = acc.getData().getValue(acc.AVERAGE.index);
        return pLat + lastDMu*latticeDensity+Math.log(sum)/volume;
    }
}
