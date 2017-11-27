/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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

    protected MCMoveStepDependent mcMove;
    protected double acceptanceTarget;
    protected long adjustInterval, maxAdjustInterval;
    protected boolean tunable = true;
    protected long lastAdjust;
    protected double adjustStep, minAdjustStep;
    protected double defaultAdjustStep = 1.05;
    protected long defaultAdjustInterval = 100;
    protected boolean noisyAdjustment = false;
    protected boolean noReset = false;

    public MCMoveStepTracker() {
        super();
        setAcceptanceTarget(0.5);
        maxAdjustInterval = Long.MAX_VALUE;
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
        double stepSize = mcMove.getStepSize();
        if (nAccept > (int) (acceptanceTarget * nTrials)) {
            if (stepSize < mcMove.getStepSizeMax()) {
                if (lastAdjust < 0) {
                    // back-and-forth
                    adjustInterval *= 2;
                    if (adjustInterval > maxAdjustInterval) adjustInterval = maxAdjustInterval;
                    adjustStep = Math.sqrt(adjustStep);
                }
                else if (lastAdjust == 5) {
                    // sixth consecutive increase.  increase adjustment step
                    adjustStep *= adjustStep;
                    if (adjustStep > 2) {
                        adjustStep = 2;
                    }
                    lastAdjust = 3;
                }
                stepSize *= adjustStep;
                stepSize = Math.min(stepSize, mcMove.getStepSizeMax());
                if (noisyAdjustment) {
                    System.out.println(mcMove+" increasing step size to "+stepSize+" (acceptance="+(double)nAccept/nTrials+")");
                }
                if (lastAdjust < 1) {
                    lastAdjust = 1;
                }
                else {
                    lastAdjust++;
                }
            }
            else if (noisyAdjustment) {
                System.out.println(mcMove+" step size: "+stepSize+" (acceptance="+(double)nAccept/nTrials+")");
            }
        } else {
            if (stepSize > mcMove.getStepSizeMin()) {
                if (lastAdjust > 0) {
                    // back-and-forth
                    adjustInterval *= 2;
                    if (adjustInterval > maxAdjustInterval) adjustInterval = maxAdjustInterval;
                    adjustStep = Math.sqrt(adjustStep);
                }
                else if (lastAdjust == -5) {
                    // sixth consecutive decrease.  increase adjustment step
                    adjustStep *= adjustStep;
                    if (adjustStep > 2) {
                        adjustStep = 2;
                    }
                    lastAdjust = -3;
                }
                stepSize /= adjustStep;
                stepSize = Math.max(stepSize, mcMove.getStepSizeMin());
                if (noisyAdjustment) {
                    System.out.println(mcMove+" decreasing step size to "+stepSize+" (acceptance="+(double)nAccept/nTrials+")");
                }
                if (lastAdjust > -1) {
                    lastAdjust = -1;
                }
                else {
                    lastAdjust--;
                }
            }
            else if (noisyAdjustment) {
                System.out.println(mcMove+" step size: "+stepSize+" (acceptance="+(double)nAccept/nTrials+")");
            }
        }
        // setStepSize will call resetAdjustStep, which is bad, so prevent it
        noReset = true;
        mcMove.setStepSize(stepSize);
        noReset = false;
        nTrials = 0;
        nAccept = 0;
        chiSum = 0;
    }

    /**
     * Resets damped step adjustment.  The adjustment interval is
     * reset to defaultAdjustInterval and the adjustment step is
     * reset to defaultAdjustStep.
     */
    public void resetAdjustStep() {
        if (noReset) {
            return;
        }
        adjustInterval = defaultAdjustInterval;
        adjustStep = defaultAdjustStep;
        lastAdjust = 0;
        chiSum = 0.0;
        nTrials = 0;
        nAccept = 0;
    }

    /**
     * The step size is adjusted every adjustInterval trials.  When the
     * direction of the change oscillates, the adjustInterval is increased,
     * but it will not be increased beyond maxAdjustInterval
     *
     * @param maxAdjustInterval the value beyond which adjustInterval will not
     *                          be increased.
     */
    public void setMaxAdjustInterval(long maxAdjustInterval) {
        this.maxAdjustInterval = maxAdjustInterval;
    }

    /**
     * When the step size is adjusted, it is adjusted by adjustStep.  When the
     * direction of the change oscillates, the adjustStep is decreased,
     * but it will not be decreased beyond minAdjustStep.
     *
     * @param minAdjustStep the value beyond which adjustStep will not be
     *                      decreased.
     */
    public void setMinAdjustStep(double minAdjustStep) {
        this.minAdjustStep = minAdjustStep;
    }

    /**
     * @return current value of the targeted rate of acceptance for this move.
     */
    public final double getAcceptanceTarget() {
        return acceptanceTarget;
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

    public final long getAdjustInterval() {
        return adjustInterval;
    }

    /**
     * Sets the interval between steps size adjustments.  This
     * also resets the acceptance averages and other step adjustment
     * parameters to their default values.
     */
    public void setAdjustInterval(long i) {
        defaultAdjustInterval = i;
        resetAdjustStep();
    }

    public double getAdjustStepSize() {
        return adjustStep;
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

    public final boolean getTunable() {
        return tunable;
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

    public boolean getNoisyAdjustment() {
        return noisyAdjustment;
    }

    public void setNoisyAdjustment(boolean isNoisy) {
        noisyAdjustment = isNoisy;
    }
}
