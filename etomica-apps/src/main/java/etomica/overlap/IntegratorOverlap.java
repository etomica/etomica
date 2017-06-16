/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.overlap;

import etomica.integrator.Integrator;
import etomica.integrator.IntegratorManagerMC;

/**
 * This integrator class manages (2) sub-integrators for an overlap
 * sampling simulation. 
 */
public class IntegratorOverlap extends IntegratorManagerMC {

    public IntegratorOverlap(Integrator[] aIntegrators) {
        super(null);
        setNumSubSteps(1000);
        for (int i=0; i<aIntegrators.length; i++) {
            addIntegrator(aIntegrators[i]);
        }
        stepFrac = new double[nIntegrators];
        totNumSubSteps = new long[nIntegrators];
        totTime = new long[nIntegrators];
        setAdjustStepFraction(true);
        setAdjustInterval(1);
        //there are no global moves
        setGlobalMoveInterval(Double.POSITIVE_INFINITY);


        // and hope nobody calls add/remove Integrators
    }
    
    /**
     * Sets the DataSource that retrieves data from both boxes and provides
     * information to the integrator about their progress, which is used to
     * decide how much to sample each system.
     */
    public void setReferenceFracSource(ReferenceFracSource fracSource) {
        this.fracSource = fracSource;
        for (int i=0; i<nIntegrators; i++) {
            if (doAdjustStepFrac) {
                stepFrac[i] = 1.0/nIntegrators;
            }
            totNumSubSteps[i] = 0;
            totTime[i] = 0;
        }
    }
    
    public void reset() {
        super.reset();
        if (fracSource != null) {
            // we just want to reset our frac-related stuff
            setReferenceFracSource(fracSource);
        }
    }

    /**
     * Sets the number of total number of integrator steps the sub-integrators
     * perform for every step of the overlap integrator.  Default value is 1000.
     */
    public void setNumSubSteps(long n) {
        numSubSteps = n;
    }
    
    /**
     * Retruns the number of total number of integrator steps the sub-integrators
     * perform for every step of the overlap integrator.  Default value is 1000.
     */
    public long getNumSubSteps() {
        return numSubSteps;
    }
    
    /**
     * Sets whether to adjust the number of relative number of steps for each
     * sub-integrator.  Default is true.
     */
    public void setAdjustStepFraction(boolean b) {
        doAdjustStepFrac = b;
    }

    /**
     * Sets whether to use "aggressive" adjustment.  With aggressive adjustment,
     * if the actual reference fraction is less than the ideal, the integrator
     * only runs the reference system. 
     */
    public void setAggressiveAdjustStepFraction(boolean b) {
        doAggressiveAdjustStepFrac = b;
    }
    
    /**
     * If true, causes the integrator to optimize how much each system runs
     * based on the wall clock time needed to run steps in each system.  If
     * false, the integrator does not consider time.
     */
    public void setDoAdjustOnTime(boolean b) {
        doAdjustStepsOnTime = b;
    }

    /**
     * Sets whether to adjust the number of relative number of steps for each
     * sub-integrator.  Default is true.
     */
    public boolean isAdjustStepFraction() {
        return doAdjustStepFrac;
    }

    /**
     * Sets how frequently the appropriate reference vs. target step fraction
     * is recalculated.
     */
    public void setAdjustInterval(long newInterval) {
        adjustInterval = newInterval;
        adjustCountdown = newInterval;
    }

    public long getAdjustInterval() {
        return adjustInterval;
    }

    // Override superclass so we can run the integrators for different lenghts
    // of time.  There are no global moves.
    public void doStepInternal() {
        int totSubSteps = 0;
        if (!doAggressiveAdjustStepFrac) {
            for (int i=0; i<nIntegrators; i++) {
                long iSubSteps = (long)(numSubSteps*stepFrac[i]+0.5);
                if (i==nIntegrators-1) {
                    iSubSteps = numSubSteps - totSubSteps;
                }
                totSubSteps += iSubSteps;
                if (doAdjustStepFrac) {
                    // if we're internally adjusting the step fractions, require at least 1%
                    iSubSteps = numSubSteps/100 + (int)(numSubSteps*(1-0.01*nIntegrators) * stepFrac[i]);
                }
                totTime[i] -= System.currentTimeMillis();
                for (int j=0; j<iSubSteps; j++) {
                    integrators[i].doStep();
                }
                totTime[i] += System.currentTimeMillis();
                totNumSubSteps[i] += iSubSteps;
            }
        }
        else {
            int iBox;
            if (stepCount < 21) {
                iBox = (int)(stepCount%2);
            }
            else {
                // stepCount starts at 1
                iBox = ((double)totNumSubSteps[0])/((stepCount-1)*numSubSteps) > stepFrac[0] ? 1 : 0;
            }
            totTime[iBox] -= System.currentTimeMillis();
            for (int j=0; j<numSubSteps; j++) {
                integrators[iBox].doStep();
            }
            totTime[iBox] += System.currentTimeMillis();
            totNumSubSteps[iBox] += numSubSteps;
        }
        if (doAdjustStepFrac && --adjustCountdown == 0) {
            double nFrac;
            if (doAdjustStepsOnTime && totTime[0] > 4*totNumSubSteps[0]/numSubSteps && totTime[1] > 4*totNumSubSteps[1]/numSubSteps) {
                double idealFrac = fracSource.getIdealRefFraction((double)totTime[0]/(totTime[0]+totTime[1]));
                double tRatio = idealFrac / (1-idealFrac);
                double nRatio = tRatio * Math.sqrt(((double)totNumSubSteps[0])*totTime[1]/totNumSubSteps[1]/totTime[0]);
                nFrac = nRatio / (1+nRatio);
            }
            else {
                nFrac = fracSource.getIdealRefFraction((double)totNumSubSteps[0]/(totNumSubSteps[0]+totNumSubSteps[1]));
            }
            if (nFrac < 0.001) nFrac = 0.001;
            else if (nFrac > 0.999) nFrac = 0.999;
            stepFrac[0] = nFrac;
            stepFrac[1] = 1.0 - nFrac;
            adjustCountdown = adjustInterval;
        }
    }

    public void setRefStepFraction(double f) {
        stepFrac[0] = f;
        stepFrac[1] = 1-f;
    }
    
    public double getIdealRefStepFraction() {
        return stepFrac[0];
    }
    
    public double getRefStepFraction() {
        return totNumSubSteps[0]/(double)(totNumSubSteps[0]+totNumSubSteps[1]);
    }
    
    public double getRefTimeFraction() {
        return totTime[0]/(double)(totTime[0]+totTime[1]);
    }
    
    private static final long serialVersionUID = 1L;
    protected final double[] stepFrac;
    protected long numSubSteps;
    protected long[] totNumSubSteps;
    protected ReferenceFracSource fracSource;
    protected boolean doAdjustStepFrac, doAggressiveAdjustStepFrac;
    protected long adjustInterval, adjustCountdown;
    protected long[] totTime;
    protected boolean doAdjustStepsOnTime;
    
    public interface ReferenceFracSource {
        public double getIdealRefFraction(double oldFrac);
    }
            
}
