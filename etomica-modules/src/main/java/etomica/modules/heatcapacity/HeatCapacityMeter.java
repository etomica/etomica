/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.heatcapacity;

import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;

/**
 * Collects everything measured by the Heat Capacity module: temperatures of the separate modes,
 * the average partition of energy among modes, the calorimetric heat capacity (from the temperature
 * rise after adding heat) and the fluctuation heat capacity (from energy fluctuations at constant T).
 * <p>
 * Calorimetry compares two averages of T.  The baseline is the average since the start of the
 * current adiabatic run, or, if heat has already been added, the average over the last half of the
 * time since that addition (the previous heat addition's "after" average).  The "after" average is
 * over the second half of the time since the latest heat addition.  Isothermal data can't
 * contribute: on switching to adiabatic, the energy freezes at one sample from the canonical
 * distribution, which moves T away from its isothermal average by an amount comparable to the
 * temperature rise being measured.
 */
public class HeatCapacityMeter implements IntegratorListener {

    /** number of integrator steps between samples */
    public static final int INTERVAL = 10;
    /** fluctuation sampling starts this much time after a reset, for equilibration */
    public static final double FLUCTUATION_EQUILIBRATION = 50;
    /** "after" and baseline averages must span at least this much time to be reported */
    public static final double MIN_WINDOW = 20;
    /** number of samples averaged for each point of the temperature histories */
    public static final int PLOT_BLOCK = 20;

    protected final HeatCapacity sim;
    protected final MeterModeEnergies modes;

    protected int stepCount;
    protected double time;
    protected double T, tTrans, tRot, tVib;

    // averages for plotting; vibration is much faster than the time between plotted points
    protected int nPlot;
    protected final double[] plotSum = new double[4];
    protected final double[] plotT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN};

    // averages since the last reset or heat addition
    protected long nAvg;
    protected double sumT, sumKETrans, sumKERot, sumKEVib, sumPEVib, sumPEInter;

    protected final BlockSeries tSeries = new BlockSeries();
    protected double windowStart;
    protected boolean heatAdded;
    protected double qLast, baselineT, baselineErr;

    protected final FluctuationBlocks fluctuations = new FluctuationBlocks();
    protected double fluctuationStart;

    public HeatCapacityMeter(HeatCapacity sim) {
        this.sim = sim;
        modes = new MeterModeEnergies(sim.box);
        T = tTrans = tRot = tVib = Double.NaN;
    }

    /** Clears all averages and the heat-capacity measurements. */
    public void reset() {
        time = 0;
        stepCount = 0;
        resetAverages();
        resetCalorimetry();
        resetFluctuations();
    }

    protected void resetAverages() {
        nAvg = 0;
        sumT = sumKETrans = sumKERot = sumKEVib = sumPEVib = sumPEInter = 0;
    }

    /** Starts a new baseline for calorimetry, forgetting about any earlier heat addition. */
    public void resetCalorimetry() {
        tSeries.reset();
        windowStart = time;
        heatAdded = false;
        qLast = baselineT = baselineErr = Double.NaN;
    }

    public void resetFluctuations() {
        fluctuations.reset();
        fluctuationStart = time + FLUCTUATION_EQUILIBRATION;
    }

    public void integratorInitialized(IntegratorEvent e) {
    }

    public void integratorStepStarted(IntegratorEvent e) {
    }

    public void integratorStepFinished(IntegratorEvent e) {
        time += sim.integrator.getTimeStep();
        if (++stepCount < INTERVAL) return;
        stepCount = 0;
        sample();
    }

    protected int nMolecules() {
        return sim.box.getMoleculeList().size();
    }

    protected int nDiatomic() {
        return sim.box.getNMolecules(sim.speciesDi);
    }

    protected int nTriatomic() {
        return sim.box.getNMolecules(sim.speciesTri);
    }

    /** Number of kinetic degrees of freedom; net momentum is zero. */
    protected int fTotal() {
        return 3 * sim.box.getLeafList().size() - 3;
    }

    protected void sample() {
        int nMol = nMolecules();
        if (nMol == 0) {
            T = tTrans = tRot = tVib = Double.NaN;
            return;
        }
        computeTemperatures();
        double ke = modes.keTrans + modes.keRot + modes.keVib;

        plotSum[0] += T;
        plotSum[1] += tTrans;
        plotSum[2] += tRot;
        plotSum[3] += tVib;
        if (++nPlot == PLOT_BLOCK) {
            for (int i = 0; i < plotSum.length; i++) {
                plotT[i] = plotSum[i] / nPlot;
                plotSum[i] = 0;
            }
            nPlot = 0;
        }

        double peVib = sim.pcBonding.getLastEnergy();
        double peInter = sim.pcPair.getLastEnergy();
        nAvg++;
        sumT += T;
        sumKETrans += modes.keTrans;
        sumKERot += modes.keRot;
        sumKEVib += modes.keVib;
        sumPEVib += peVib;
        sumPEInter += peInter;

        if (sim.integrator.isIsothermal()) {
            if (time >= fluctuationStart) fluctuations.add(ke + peVib + peInter);
        } else {
            tSeries.add(time, T);
        }
    }

    private void computeTemperatures() {
        int nMol = nMolecules();
        modes.compute();
        int fTrans = 3 * nMol - 3;
        int fRot = 2 * nDiatomic() + 3 * nTriatomic();
        int fVib = nDiatomic() + 3 * nTriatomic();
        T = 2 * (modes.keTrans + modes.keRot + modes.keVib) / fTotal();
        tTrans = fTrans > 0 ? 2 * modes.keTrans / fTrans : Double.NaN;
        tRot = fRot > 0 ? 2 * modes.keRot / fRot : Double.NaN;
        tVib = fVib > 0 ? 2 * modes.keVib / fVib : Double.NaN;
    }

    /**
     * Starts the temperature histories over, from the current velocities.  The overall value
     * is always finite (even without molecules), because a plot can't draw without any data.
     */
    public void restartPlot() {
        nPlot = 0;
        java.util.Arrays.fill(plotSum, 0);
        double T0 = sim.integrator.getTemperature();
        if (nMolecules() > 1 && sim.integrator.isInitialized()) computeTemperatures();
        else T = tTrans = tRot = tVib = Double.NaN;
        plotT[0] = Double.isNaN(T) ? T0 : T;
        plotT[1] = tTrans;
        plotT[2] = tRot;
        plotT[3] = tVib;
    }

    /**
     * Adds heat to the translational motion of the molecules.
     *
     * @param amount heat per molecule, in units of the current kT
     * @return false if no heat could be added
     */
    public boolean addHeat(double amount) {
        int nMol = nMolecules();
        if (nMol < 2 || amount <= 0) return false;
        double[] baseline = tSeries.average(windowStart + (heatAdded ? 0.5 * (time - windowStart) : 0));
        double spanned = time - windowStart;
        double tRef = !Double.isNaN(baseline[0]) ? baseline[0] : (nAvg > 0 ? T : sim.integrator.getTemperature());
        double q = amount * tRef * nMol;
        if (!sim.addHeatToTranslation(q)) return false;
        boolean goodBaseline = !Double.isNaN(baseline[0]) && (heatAdded ? 0.5 * spanned : spanned) >= MIN_WINDOW;
        baselineT = goodBaseline ? baseline[0] : Double.NaN;
        baselineErr = goodBaseline ? baseline[1] : Double.NaN;
        qLast = q;
        heatAdded = true;
        windowStart = time;
        tSeries.reset();
        resetAverages();
        return true;
    }

    public double getTemperature() {
        return T;
    }

    public double getTTranslation() {
        return tTrans;
    }

    public double getTRotation() {
        return tRot;
    }

    public double getTVibration() {
        return tVib;
    }

    /**
     * @param i 0 for overall, 1 for translation, 2 for rotation, 3 for vibration
     * @return the temperature averaged over the latest block of PLOT_BLOCK samples
     */
    public double getPlotT(int i) {
        return plotT[i];
    }

    /** Average kinetic temperature since the last reset or heat addition */
    public double getAverageT() {
        return nAvg == 0 ? Double.NaN : sumT / nAvg;
    }

    public double getElapsedTime() {
        return time;
    }

    public boolean hasBaseline() {
        return !Double.isNaN(baselineT);
    }

    public boolean isHeatAdded() {
        return heatAdded;
    }

    /** Time since the most recent heat addition. */
    public double getTimeSinceHeat() {
        return heatAdded ? time - windowStart : Double.NaN;
    }

    /**
     * @return average per molecule of translational KE, rotational KE, vibrational KE, vibrational
     * PE and intermolecular PE, in units of the average kT
     */
    public double[] getPartition() {
        double[] p = new double[5];
        int nMol = nMolecules();
        if (nAvg == 0 || nMol == 0) {
            java.util.Arrays.fill(p, Double.NaN);
            return p;
        }
        double norm = 1.0 / (nMol * sumT);
        p[0] = sumKETrans * norm;
        p[1] = sumKERot * norm;
        p[2] = sumKEVib * norm;
        p[3] = sumPEVib * norm;
        p[4] = sumPEInter * norm;
        return p;
    }

    /** @return the equipartition values corresponding to getPartition, for the current composition */
    public double[] getEquipartition() {
        double[] p = new double[5];
        int nMol = nMolecules();
        if (nMol == 0) {
            java.util.Arrays.fill(p, Double.NaN);
            return p;
        }
        double nd = nDiatomic(), nt = nTriatomic();
        p[0] = 1.5;
        p[1] = (nd + 1.5 * nt) / nMol;
        p[2] = (0.5 * nd + 1.5 * nt) / nMol;
        p[3] = (0.5 * nd + 1.5 * nt) / nMol;
        p[4] = 0;
        return p;
    }

    /** @return the mole-fraction-weighted classical Cv/k per molecule */
    public double getPredictedCv() {
        double s = 0;
        for (double v : getEquipartition()) s += v;
        return s;
    }

    /** @return {Cv/k per molecule, its error} from the latest heat addition; NaN if not available */
    public double[] getCalorimetricCv() {
        double[] r = {Double.NaN, Double.NaN};
        if (!heatAdded || !hasBaseline()) return r;
        double dt = time - windowStart;
        if (0.5 * dt < MIN_WINDOW) return r;
        double[] after = tSeries.average(windowStart + 0.5 * dt);
        if (Double.isNaN(after[0])) return r;
        double dT = after[0] - baselineT;
        r[0] = qLast / (nMolecules() * dT);
        r[1] = Math.abs(r[0]) * Math.sqrt(after[1] * after[1] + baselineErr * baselineErr) / Math.abs(dT);
        return r;
    }

    /** @return {Cv/k per molecule, its error} from energy fluctuations; NaN if not available */
    public double[] getFluctuationCv() {
        double T0 = sim.integrator.getTemperature();
        double[] var = fluctuations.variance();
        double norm = 1.0 / (T0 * T0 * nMolecules());
        return new double[]{var[0] * norm, var[1] * norm};
    }

    /**
     * Series of block averages with fixed storage.  Blocks are merged in pairs when storage fills,
     * so the blocks always have equal size and span the whole series.
     */
    static class BlockSeries {
        private static final int CAPACITY = 1024;
        private static final int N_BINS = 10;
        private final double[] blockTime = new double[CAPACITY];
        private final double[] blockValue = new double[CAPACITY];
        private int nBlocks, blockSize = 1, nInBlock;
        private double sum;

        void reset() {
            nBlocks = 0;
            blockSize = 1;
            nInBlock = 0;
            sum = 0;
        }

        void add(double t, double v) {
            sum += v;
            if (++nInBlock < blockSize) return;
            blockTime[nBlocks] = t;
            blockValue[nBlocks] = sum / nInBlock;
            nBlocks++;
            sum = 0;
            nInBlock = 0;
            if (nBlocks == CAPACITY) {
                for (int i = 0; i < CAPACITY / 2; i++) {
                    blockTime[i] = blockTime[2 * i + 1];
                    blockValue[i] = 0.5 * (blockValue[2 * i] + blockValue[2 * i + 1]);
                }
                nBlocks = CAPACITY / 2;
                blockSize *= 2;
            }
        }

        /**
         * @return {average, error} of the blocks finishing at or after the given time.  The error
         * comes from rebinning into 10 bins and is NaN if there are not enough blocks for that.
         */
        double[] average(double tStart) {
            int first = 0;
            while (first < nBlocks && blockTime[first] < tStart) first++;
            int n = nBlocks - first;
            if (n == 0) return new double[]{Double.NaN, Double.NaN};
            double s = 0;
            for (int i = first; i < nBlocks; i++) s += blockValue[i];
            double avg = s / n;
            if (n < N_BINS) return new double[]{avg, Double.NaN};
            double s2 = 0;
            for (int b = 0; b < N_BINS; b++) {
                int i0 = first + (int) ((long) n * b / N_BINS);
                int i1 = first + (int) ((long) n * (b + 1) / N_BINS);
                double binSum = 0;
                for (int i = i0; i < i1; i++) binSum += blockValue[i];
                double d = binSum / (i1 - i0) - avg;
                s2 += d * d;
            }
            return new double[]{avg, Math.sqrt(s2 / (N_BINS * (N_BINS - 1)))};
        }
    }

    /**
     * Variance of a fluctuating quantity with a block-averaged error.  Block length doubles as
     * needed, so that for a long run the blocks become much longer than the correlation time.
     */
    static class FluctuationBlocks {
        private static final int CAPACITY = 64;
        private static final int MIN_BLOCKS = CAPACITY / 2;
        private final double[] blockX = new double[CAPACITY];
        private final double[] blockX2 = new double[CAPACITY];
        private int nBlocks, blockSize = 1, nInBlock;
        private double sumX, sumX2, x0;
        private boolean started;

        void reset() {
            nBlocks = 0;
            blockSize = 1;
            nInBlock = 0;
            sumX = sumX2 = 0;
            started = false;
        }

        void add(double v) {
            if (!started) {
                // measure relative to first value, to avoid taking a small difference of large numbers
                x0 = v;
                started = true;
            }
            double x = v - x0;
            sumX += x;
            sumX2 += x * x;
            if (++nInBlock < blockSize) return;
            blockX[nBlocks] = sumX / nInBlock;
            blockX2[nBlocks] = sumX2 / nInBlock;
            nBlocks++;
            sumX = sumX2 = 0;
            nInBlock = 0;
            if (nBlocks == CAPACITY) {
                for (int i = 0; i < CAPACITY / 2; i++) {
                    blockX[i] = 0.5 * (blockX[2 * i] + blockX[2 * i + 1]);
                    blockX2[i] = 0.5 * (blockX2[2 * i] + blockX2[2 * i + 1]);
                }
                nBlocks = CAPACITY / 2;
                blockSize *= 2;
            }
        }

        /** @return {variance, error}; both are NaN until the blocks have been merged at least 5 times. */
        double[] variance() {
            if (blockSize < 32 || nBlocks < MIN_BLOCKS) return new double[]{Double.NaN, Double.NaN};
            double avgX = 0, avgX2 = 0;
            for (int i = 0; i < nBlocks; i++) {
                avgX += blockX[i];
                avgX2 += blockX2[i];
            }
            avgX /= nBlocks;
            avgX2 /= nBlocks;
            // first-order propagation:  each block contributes X2 - 2 <X> X
            double avgG = avgX2 - 2 * avgX * avgX;
            double s2 = 0;
            for (int i = 0; i < nBlocks; i++) {
                double d = blockX2[i] - 2 * avgX * blockX[i] - avgG;
                s2 += d * d;
            }
            return new double[]{avgX2 - avgX * avgX, Math.sqrt(s2 / (nBlocks * (nBlocks - 1)))};
        }
    }
}
