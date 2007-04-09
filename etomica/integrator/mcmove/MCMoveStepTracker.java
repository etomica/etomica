package etomica.integrator.mcmove;

/**
 * This Class tracks acceptance statistics for an MCMove and also manages the
 * move's step size.  Given an acceptance target, it can adjust the 
 * step size (subject to min/max constraints) to achieve the desired acceptance
 * ratio.
 *
 * @author Andrew Schultz
 */
public class MCMoveStepTracker extends MCMoveTracker {

    public MCMoveStepTracker() {
        super();
        setAcceptanceTarget(0.5);
    }
    
    public void setMCMove(MCMoveStepDependent newMCMove) {
        mcMove = newMCMove;
    }
    
    public void updateCounts(boolean moveWasAccepted, double chi) {
        super.updateCounts(moveWasAccepted, chi);
        if (tunable && nTrials > adjustInterval) {
            adjustStepSize();
        }
    }

    protected void adjustStepSize() {
        if (nTrials == 0) {
            return;
        }
        nTrialsSum += nTrials;
        nAcceptSum += nAccept;
        double stepSize = mcMove.getStepSize();
        if (nAccept > (int) (acceptanceTarget * nTrials)) {
            if (stepSize < mcMove.getStepSizeMax()) {
                if (lastAdjust == -1) {
                    // back-and-forth
                    adjustInterval *= 2;
                    adjustStep *= 0.5;
                }
                stepSize *= 1.0+adjustStep;
                if (noisyAdjustment) {
                    System.out.println(mcMove.getClass()+" increasing step size to "+stepSize+" (acceptance="+(double)nAccept/nTrials+")");
                }
                lastAdjust = 1;
            }
        } else {
            if (stepSize > mcMove.getStepSizeMin()) {
                if (lastAdjust == 1) {
                    // back-and-forth
                    adjustInterval *= 2;
                    adjustStep *= 0.5;
                }
                stepSize *= 1.0-adjustStep;
                if (noisyAdjustment) {
                    System.out.println(mcMove.getClass()+" decreasing step size to "+stepSize+" (acceptance="+(double)nAccept/nTrials+")");
                }
                lastAdjust = -1;
            }
        }
        stepSize = Math.min(stepSize, mcMove.getStepSizeMax());
        stepSize = Math.max(stepSize, mcMove.getStepSizeMin());
        // setStepSize will call resetAdjustStep, which is bad, so prevent it
        noReset = true;
        mcMove.setStepSize(stepSize);
        noReset = false;
        nTrials = 0;
        nAccept = 0;
    }

    /**
     * Sets the desired rate of acceptance (as a fraction between 0 and 1
     * inclusive) of trials of this move. IllegalArgumentException is thrown if
     * given value is outside of acceptable range.  The average acceptance 
     * ratio and probability and the step size adjustment parameters and 
     * convergence process are all reset.
     */
    public final void setAcceptanceTarget(double target) {
        if (target < 0.0 || target > 1.0) {
            throw new IllegalArgumentException(
                    "Acceptance target should be a number between zero and unity");
        }
        acceptanceTarget = target;
        resetAdjustStep();
    }

    /**
     * Resets damped step adjustment.  The adjustment interval is
     * reset to defaultAdjustInterval and the adjustmentment step is
     * reset to defaultAdjustStep.
     */
    public void resetAdjustStep() {
        if (noReset) {
            return;
        }
        adjustInterval = defaultAdjustInterval;
        adjustStep = defaultAdjustStep;
        lastAdjust = 0;
        nTrialsSum = 0;
        nAcceptSum = 0;
        chiSum = 0.0;
        nTrials = 0;
        nAccept = 0;
    }
    
    /**
     * @return current value of the targeted rate of acceptance for this move.
     */
    public final double getAcceptanceTarget() {
        return acceptanceTarget;
    }

    /**
     * Sets the interval between steps size adjustments.  This
     * also resets the acceptance averages and other step adjustment 
     * parameters to their default values.
     */
    public void setAdjustInterval(int i) {
        defaultAdjustInterval = i;
        resetAdjustStep();
    }

    public final int getAdjustInterval() {
        return adjustInterval;
    }

    /**
     * Sets the step adjustment size.  The step size is multiplied by 
     * (1+s) when the acceptance ratio is too high and (1-s) when the step 
     * size is too low.  This also resets the acceptance averages and other 
     * step adjustment parameters to their default values.
     */
    public void setAdjustStepSize(double s) {
        defaultAdjustStep = s;
        resetAdjustStep();
    }
    
    public double getAdjustStepSize() {
        return adjustStep;
    }
    
    /**
     * Sets a flag to indicate whether tuning of the move is to be performed
     * Tuning adjust the step size or other property of the move with the aim of
     * achieving a trial acceptance rate that is near to some target value. Some
     * moves (e.g., simple insertion trial) are inherently untunable, and will
     * have false for this property.
     */
    public final void setTunable(boolean b) {
        tunable = b;
    }

    public final boolean getTunable() {
        return tunable;
    }
    
    public void setNoisyAdjustment(boolean isNoisy) {
        noisyAdjustment = isNoisy;
    }
    
    public boolean getNoisyAdjustment() {
        return noisyAdjustment;
    }

    private static final long serialVersionUID = 1L;
    protected MCMoveStepDependent mcMove;
    protected double acceptanceTarget;
    protected int nTrialsSum, nAcceptSum, adjustInterval;
    protected boolean tunable = true;
    protected int lastAdjust;
    protected double adjustStep;
    protected double defaultAdjustStep = 0.05;
    protected int defaultAdjustInterval = 100;
    protected boolean noisyAdjustment = false;
    protected boolean noReset = false;
}
