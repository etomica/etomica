/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataSourceScalar;
import etomica.integrator.mcmove.MCMoveOverlapListener;
import etomica.units.Pressure;

/**
 * Computes pressure from the grand canonical route based on measured free
 * energy differences and simulation parameters.
 * 
 * @author Andrew Schultz
 */
public class DataSourceAvgPressure2 extends DataSourceScalar {
    protected final MCMoveOverlapListener mcMoveOverlapMeter;
    protected double bmu;
    protected double bALattice;
    protected double latticeDensity;
    protected double volume;

    public DataSourceAvgPressure2(MCMoveOverlapListener mcMoveOverlapMeter, double betaMu, double betaALattice, double latticeDensity, double volume) {
        super("pressure", Pressure.DIMENSION);
        this.mcMoveOverlapMeter = mcMoveOverlapMeter;
        this.bmu = betaMu;
        this.bALattice = betaALattice;
        this.latticeDensity = latticeDensity;
        this.volume = volume;
    }

    public void setBetaMu(double newBetaMu) {
        bmu = newBetaMu;
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
        return (bmu-bALattice)*latticeDensity+Math.log(sum)/volume;
    }
}