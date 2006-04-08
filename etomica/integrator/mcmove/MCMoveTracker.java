package etomica.integrator.mcmove;

import java.io.Serializable;

/**
 * This Class is responsible for tracking acceptance statistics for an MCMove.
 * The acceptance ratio and probability can be retrieved.  This class has no 
 * effect on the MCMove.
 *
 * @author Andrew Schultz
 */
public class MCMoveTracker implements Serializable {

    /**
     * Updates statistics regarding the acceptance rate of this move.  This 
     * method should only be called by the Integrator after a move has been
     * accepted or rejected.
     * 
     * @param moveWasAccepted
     *            indicates whether the most recently attempted move was
     *            accepted
     * @param chi
     *            indicates the probability with which the MC move was
     *            accepted.
     */
    public void updateCounts(boolean moveWasAccepted, double chi) {
        nTrials++;
        if (moveWasAccepted)
            nAccept++;
        chiSum += chi;
    }

    /**
     * Returns the fraction of accepted trials of the move associated with 
     * this tracker.
     */
    public double acceptanceRatio() {
        return (nTrials > 0) ? (double) nAccept / (double) nTrials
                        : Double.NaN;
    }


    /**
     * Returns the average probability of accepting the move associated with
     * this tracker.
     */
    public double acceptanceProbability() {
        return ((nTrials > 0) ? chiSum / nTrials : Double.NaN);
    }
    
    protected int nTrials, nAccept;
    protected double chiSum;
}
