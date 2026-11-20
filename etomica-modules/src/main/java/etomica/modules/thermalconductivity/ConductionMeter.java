/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.thermalconductivity;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;

import java.util.ArrayList;
import java.util.List;

/**
 * Collects the quantities needed to measure the thermal conductivity: binned kinetic
 * temperature profile, heat fluxes delivered/removed by the walls, the linear fit
 * to the interior of the profile, and k = J/(-dT/dx) with a block-average error.
 * Runs as an integrator listener (every step).  All accumulators restart on reset().
 */
public class ConductionMeter implements IntegratorListener {

    public static final int NBINS = 20;
    /** number of bins excluded at each wall from the linear fit */
    public static final int N_EXCLUDE = 2;
    /** minimum samples in a bin for its temperature to be reported */
    protected static final double MIN_SAMPLES = 20;

    protected final ThermalConductivity sim;
    protected final double[] sumKE = new double[NBINS], sumN = new double[NBINS];
    protected final double[] blkKE = new double[NBINS], blkN = new double[NBINS];
    protected double blockTime = 20;
    protected double tStart, tBlockStart, lastT;
    protected double heatHotBlockStart, heatColdBlockStart;

    protected final List<double[]> blocks = new ArrayList<>(); // {tEnd, Jhot, Jcold, dT/dx}

    public ConductionMeter(ThermalConductivity sim) {
        this.sim = sim;
    }

    public synchronized void reset() {
        double t = sim.integrator.getCurrentTime();
        resetAt(t);
    }

    private void resetAt(double t) {
        java.util.Arrays.fill(sumKE, 0);
        java.util.Arrays.fill(sumN, 0);
        java.util.Arrays.fill(blkKE, 0);
        java.util.Arrays.fill(blkN, 0);
        blocks.clear();
        sim.wall.resetHeat();
        heatHotBlockStart = heatColdBlockStart = 0;
        tStart = tBlockStart = lastT = t;
    }

    public void setBlockTime(double blockTime) {
        this.blockTime = blockTime;
    }

    @Override
    public synchronized void integratorStepFinished(IntegratorEvent e) {
        double t = sim.integrator.getCurrentTime();
        if (t < lastT) {
            // integrator time was restarted
            resetAt(t);
        }
        lastT = t;
        double lx = sim.getLx();
        double m = sim.species.getLeafType().getMass();
        IAtomList atoms = sim.box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            IAtomKinetic a = (IAtomKinetic) atoms.get(i);
            double x = a.getPosition().getX(0) + 0.5 * lx;
            int bin = (int) (x / lx * NBINS);
            if (bin < 0) bin = 0;
            else if (bin >= NBINS) bin = NBINS - 1;
            double ke = 0.5 * m * a.getVelocity().squared();
            sumKE[bin] += ke;
            sumN[bin]++;
            blkKE[bin] += ke;
            blkN[bin]++;
        }
        if (t - tBlockStart >= blockTime) {
            closeBlock(t);
        }
    }

    protected void closeBlock(double t) {
        double dt = t - tBlockStart;
        double ly = sim.getLy();
        double jHot = (sim.wall.heatHot - heatHotBlockStart) / (dt * ly);
        double jCold = -(sim.wall.heatCold - heatColdBlockStart) / (dt * ly);
        double[] tProfile = profile(blkKE, blkN);
        double[] fit = fitProfile(tProfile);
        double slope = fit == null ? Double.NaN : fit[0];
        blocks.add(new double[]{t - tStart, jHot, jCold, slope});
        heatHotBlockStart = sim.wall.heatHot;
        heatColdBlockStart = sim.wall.heatCold;
        java.util.Arrays.fill(blkKE, 0);
        java.util.Arrays.fill(blkN, 0);
        tBlockStart = t;
    }

    protected double[] profile(double[] ke, double[] n) {
        double[] T = new double[NBINS];
        for (int i = 0; i < NBINS; i++) {
            // T = <KE per disk> / (D/2 kB), D = 2
            T[i] = n[i] >= MIN_SAMPLES ? ke[i] / n[i] / (sim.getSpace().D() / 2.0) : Double.NaN;
        }
        return T;
    }

    /** x position (0..Lx) of the center of each bin */
    public double[] binCenters() {
        double[] x = new double[NBINS];
        double dx = sim.getLx() / NBINS;
        for (int i = 0; i < NBINS; i++) x[i] = (i + 0.5) * dx;
        return x;
    }

    /** Temperature profile accumulated since last reset (NaN for bins with too few samples). */
    public synchronized double[] getProfile() {
        return profile(sumKE, sumN);
    }

    /**
     * Least-squares line through the interior bins of the given profile.
     *
     * @return {slope, intercept (T at x=0), mean T of the bins used, number of bins used}
     * or null if fewer than 3 bins are usable.
     */
    protected double[] fitProfile(double[] T) {
        double[] x = binCenters();
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int n = 0;
        for (int i = N_EXCLUDE; i < NBINS - N_EXCLUDE; i++) {
            if (Double.isNaN(T[i])) continue;
            n++;
            sx += x[i];
            sy += T[i];
            sxx += x[i] * x[i];
            sxy += x[i] * T[i];
        }
        if (n < 3) return null;
        double denom = n * sxx - sx * sx;
        double slope = (n * sxy - sx * sy) / denom;
        double intercept = (sy - slope * sx) / n;
        return new double[]{slope, intercept, sy / n, n};
    }

    /** Fit of the cumulative profile: {slope, intercept, mean T, n} or null */
    public synchronized double[] getFit() {
        return fitProfile(getProfile());
    }

    public synchronized double getElapsedTime() {
        return lastT - tStart;
    }

    public synchronized double getJHot() {
        double dt = lastT - tStart;
        return dt > 0 ? sim.wall.heatHot / (dt * sim.getLy()) : Double.NaN;
    }

    public synchronized double getJCold() {
        double dt = lastT - tStart;
        return dt > 0 ? -sim.wall.heatCold / (dt * sim.getLy()) : Double.NaN;
    }

    public synchronized double getJMean() {
        return 0.5 * (getJHot() + getJCold());
    }

    /** @return {k, error, T at midpoint} (k = NaN if not available; error = NaN with fewer than 3 blocks) */
    public synchronized double[] getConductivity() {
        double[] fit = getFit();
        double[] result = {Double.NaN, Double.NaN, Double.NaN};
        if (fit == null || fit[0] >= 0 || Double.isNaN(getJMean())) return result;
        result[0] = getJMean() / (-fit[0]);
        result[2] = fit[2];
        // block-average error using the linearization dk/k = dJ/J - ds/s, which stays
        // well behaved when an individual block has a slope near zero
        double jBar = getJMean();
        double sBar = fit[0];
        double s = 0, s2 = 0;
        int nb = 0;
        for (double[] b : blocks) {
            if (Double.isNaN(b[3])) continue;
            double f = (0.5 * (b[1] + b[2]) - jBar) / jBar - (b[3] - sBar) / sBar;
            nb++;
            s += f;
            s2 += f * f;
        }
        if (nb >= 3) {
            double var = (s2 - s * s / nb) / (nb - 1);
            result[1] = Math.abs(result[0]) * Math.sqrt(Math.max(var, 0) / nb);
        }
        return result;
    }

    /** Block history as {tEnd, Jhot, Jcold, dT/dx} rows (copied). */
    public synchronized double[][] getBlocks() {
        double[][] out = new double[blocks.size()][];
        for (int i = 0; i < out.length; i++) out[i] = blocks.get(i).clone();
        return out;
    }

    public void printSummary() {
        double[] T = getProfile();
        double[] x = binCenters();
        for (int i = 0; i < NBINS; i++) {
            System.out.printf("%6.2f %8.4f%n", x[i], T[i]);
        }
        double[] kk = getConductivity();
        double[] fit = getFit();
        System.out.printf("phi %5.3f N %d  Jhot %9.5f Jcold %9.5f  dT/dx %9.5f  k %8.4f +/- %8.4f  k/sqrtT %8.4f (dilute CE %8.4f)%n",
                sim.getPhi(), sim.box.getLeafList().size(), getJHot(), getJCold(),
                fit == null ? Double.NaN : fit[0], kk[0], kk[1], kk[0] / Math.sqrt(kk[2]), 2 / Math.sqrt(Math.PI));
    }
}
