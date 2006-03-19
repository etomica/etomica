package etomica.integrator;

import etomica.atom.iterator.AtomIterator;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;

/**
 * Parent class for all elementary Monte Carlo move classes, as used by
 * IntegratorMC class. Includes actions needed to complete trial move, and
 * reports numbers needed to determine if move is accepted. <br>
 * Subclasses should set value of nominal frequency and perParticleFrequency
 * flag if the move is to be performed at a rate different from the defaults.
 * Default nominalFrequency is 100, and perParticleFrequency is false.
 * 
 * @author David Kofke
 */

public abstract class MCMove implements java.io.Serializable {

    /**
     * @param potentialMaster the potential master that move can use to calculate energy
     * @param nPhases the number of phases on which the move acts.  This is used at
     * construction to size the (final) phases array and cannot be changed.
     */
	public MCMove(PotentialMaster potentialMaster, int nPhases) {
        this.potential = potentialMaster;
        setAcceptanceTarget(0.5);
        nominalFrequency = 100;
        perParticleFrequency = false;
    }

	/**
	 * Updates statistics regarding the acceptance rate of this move.
	 * 
	 * @param moveWasAccepted
	 *            indicates whether the most recently attempted move was
	 *            accepted
	 * @param doAdjustStepSize
	 *            indicates whether move is permitted to perform any adjustments
	 *            of its step size based on the acceptance history; typically
	 *            this will be false if the integrator is in production mode,
	 *            and true if it is still equilibrating. Adjustment will be
	 *            performed only if doAdjustStepSize is true and current number
	 *            of trials since last adjustment exceeds adjustInterval.
	 */
	public void updateCounts(boolean moveWasAccepted, double chi, boolean doAdjustStepSize) {
		nTrials++;
		if (moveWasAccepted)
			nAccept++;
        chiSum += chi;
		if (doAdjustStepSize && nTrials > adjustInterval) {
			adjustStepSize();
		}
	}

	/**
	 * Method to perform trial move. Returns false if the trial could not be
	 * attempted, for example if there were no molecules in the phase and the
	 * trial is designed to displace an atom; returns true otherwise.
	 */
	public abstract boolean doTrial();

	/**
     * Returns the temperature-independent part of the acceptance probability.  
     * The actual acceptance probility is calculated as max(1,A*exp(-B/T))
     * where T is the IntegratorPhase's temperature.
	 */
	public abstract double getA();

	/**
	 * Returns the temperature-dependent part of the acceptance probability.  
     * The actual acceptance probility is calculated as max(1,A*exp(-B/T))
     * where T is the IntegratorPhase's temperature.
	 */
	public abstract double getB();

	/**
	 * Method called by IntegratorMC in the event that the most recent trial is
	 * accepted.
	 */
	public abstract void acceptNotify();

	/**
	 * Method called by IntegratorMC in the event that the most recent trial
	 * move is rejected. This method should cause the system to be restored to
	 * the condition before the most recent call to doTrial.
	 */
	public abstract void rejectNotify();

	/**
	 * Returns an iterator that yields the atoms that were affected by the trial
	 * move the last time doTrial was invoked (regardless of whether the move
	 * was accepted). This information usually is not needed, but it is
	 * available in cases where required by objects that perform some activity
	 * while the move is in progress, or need to update after the move is
	 * completed. Such objects can receive notification of the move's completion
	 * by registering with the IntegratorMC as MCMoveEventListeners.  If the 
     * move caused an atom to be removed from the phase, it will not be returned
     * because that Atom is not considered an affected Atom.
	 */
	public abstract AtomIterator affectedAtoms(Phase phase);

	/**
	 * Returns the change in the energy of the given phase that results from the
	 * trial move. Should be called only after lnProbabilityRatio(); returns
	 * Double.NaN if invoked between calls to doTrial and lnProbabilityRatio.
	 */
	public abstract double energyChange(Phase phase);

    /**
     * Sets the phase(s) on which this move acts.  The number of phases is set
     * at construction, and the size of the given array must equal this value.
     * The phase(s) itself can be changed via this method, if desired.  Most moves
     * act on only one phase, but some (such as a volume exchange move) involve
     * two or more phases.
     * @param p
     */
	public void setPhase(Phase p) {
        phase = p;
	}

    /**
     * @return the phase(s) on which this move acts.
     */
	public Phase getPhase() {
		return phase;
	}

    protected void adjustStepSize() {
        if (nTrials == 0) {
            return;
        }
        nTrialsSum += nTrials;
        nAcceptSum += nAccept;
        if (nTrialsSum != 0) {
            acceptanceRatio = (double) nAcceptSum / (double) nTrialsSum;
            acceptanceProbability = chiSum / nTrialsSum;
        }
        if (!tunable) {
            nTrials = nAccept = 0;
            return;
        }
        if (nAccept > (int) (acceptanceTarget * nTrials)) {
            if (stepSize < stepSizeMax) {
                if (lastAdjust == -1) {
                    // back-and-forth
                    adjustInterval *= 2;
                    adjustStep *= 0.5;
                }
                stepSize *= 1.0+adjustStep;
                if (noisyAdjustment) {
                    System.out.println(phase+" "+this.getClass()+" increasing step size to "+stepSize+" (acceptance="+(double)nAccept/nTrials+")");
                }
                lastAdjust = 1;
            }
        } else {
            if (stepSize > stepSizeMin) {
                if (lastAdjust == 1) {
                    // back-and-forth
                    adjustInterval *= 2;
                    adjustStep *= 0.5;
                }
                stepSize *= 1.0-adjustStep;
                if (noisyAdjustment) {
                    System.out.println(phase+" "+this.getClass()+" decreasing step size to "+stepSize+" (acceptance="+(double)nAccept/nTrials+")");
                }
                lastAdjust = -1;
            }
        }
        stepSize = Math.min(stepSize, stepSizeMax);
        stepSize = Math.max(stepSize, stepSizeMin);
        nTrials = 0;
        nAccept = 0;
    }

    /**
     * Sets whether to print out information about step size whenever it
     * is adjusted.
     */
    public void setNoisyAdjustment(boolean noisy) {
        noisyAdjustment = noisy;
    }
    
    /**
     * Resets damped step adjustment.  The adjustment interval is
     * reset to defaultAdjustInterval and the adjustmentment step is
     * reset to defaultAdjustStep.
     */
    public void resetAdjustStep() {
        adjustInterval = defaultAdjustInterval;
        adjustStep = defaultAdjustStep;
        lastAdjust = 0;
        nTrialsSum = 0;
        nAcceptSum = 0;
        chiSum = 0.0;
        nTrials = 0;
        nAccept = 0;
        acceptanceRatio = Double.NaN;
        acceptanceProbability = Double.NaN;
    }
    
	/**
	 * Returns a nominal, unnormalized frequency for performing this move,
	 * relative to the other moves that have been added to the integrator. Each
	 * move is performed (on average) an amount in proportion to this frequency.
	 * Moves having the same frequency are performed with equal likelihood.
	 * Default value is 100, but this may be overridden by the subclass. The
	 * nominal frequency cannot be changed, as it is used only to advise the
	 * integrator when the move is added to it. Adjustment of the move frequency
	 * (if the nominal value is not desired) can be accomplished via the
	 * setFrequency method of IntegratorMC.
	 */
	public final int getNominalFrequency() {
		return nominalFrequency;
	}

	/**
	 * Indicates whether this move should nominally be performed at a frequency
	 * proportional to the number of molecules in the phase.
	 * 
	 * @see #nominalFrequency
	 */
	public final boolean isNominallyPerParticleFrequency() {
		return perParticleFrequency;
	}

	/**
	 * Fraction of time trials of this type were accepted since acceptanceTarget
	 * was set.
	 */
	public double acceptanceRatio() {
		return (acceptanceRatio >= 0.0) ? acceptanceRatio
				: (nTrials > 0) ? (double) nAccept / (double) nTrials
						: Double.NaN;
	}


    /**
     * Fraction of time trials of this type were accepted since acceptanceTarget
     * was set.
     */
    public double acceptanceProbability() {
        return (acceptanceProbability >= 0.0) ? acceptanceProbability :
                ((nTrials > 0) ? chiSum / nTrials : Double.NaN);
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
	 * @return current value of the targeted rate of acceptance for this move.
	 */
	public final double getAcceptanceTarget() {
		return acceptanceTarget;
	}

    /**
     * Sets the step size of the move.  The average acceptance ratio and
     * probability are reset along with the step size adjustment parameters.  
     */
	public void setStepSize(double step) {
		stepSize = step;
        resetAdjustStep();
	}

    /**
     * returns the current step size.
     */
	public double getStepSize() {
		return stepSize;
	}

	public void setStepSizeMax(double step) {
		stepSizeMax = step;
	}

	public double getStepSizeMax() {
		return stepSizeMax;
	}

	public void setStepSizeMin(double step) {
		stepSizeMin = step;
	}

	public double getStepSizeMin() {
		return stepSizeMin;
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

	/**
	 * Accessor method of the name of this object
	 * 
	 * @return The given name
	 */
	public final String getName() {
		return name;
	}

	/**
	 * Method to set the name of this object
	 * 
	 * @param name
	 *            The name string to be associated with this object
	 */
	public final void setName(String name) {
		this.name = name;
	}

	/**
	 * Overrides the Object class toString method to have it return the output
	 * of getName
	 * 
	 * @return The name given to the object
	 */
	public String toString() {
		return getName();
	}

	/**
	 * Value giving nominal frequency for performing this move. Default is 100,
	 * but may be given a different value by subclasses.
	 */
	protected int nominalFrequency;

	/**
	 * Flag indicating whether nominal frequency is interpreted as a
	 * perParticleFrequency, or as a full frequency. Default is false, but may
	 * be given a different value by subclasses.
	 */
	protected boolean perParticleFrequency;
	private double acceptanceRatio, acceptanceTarget, acceptanceProbability;
	protected double stepSize, stepSizeMax, stepSizeMin;
	private int nTrials, nAccept, nTrialsSum, nAcceptSum, adjustInterval;
    private double chiSum;
	protected boolean tunable = true;
	protected Phase phase;
	private String name;
	protected final PotentialMaster potential;
    private int lastAdjust;
    private double adjustStep;
    private double defaultAdjustStep = 0.05;
    private int defaultAdjustInterval = 100;
    private boolean noisyAdjustment = false;

}