package etomica;

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
		nTrials = 0;
		nAccept = 0;
		setAcceptanceTarget(0.5);
		nominalFrequency = 100;
		perParticleFrequency = false;
		setAdjustInterval(100);
        phases = new Phase[nPhases];
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
	public void updateCounts(boolean moveWasAccepted, boolean doAdjustStepSize) {
		nTrials++;
		if (moveWasAccepted)
			nAccept++;
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
	 * Returns log of the ratio of the trial probabilities, ln(Tji/Tij) for the
	 * states encountered before (i) and after (j) the most recent call to
	 * doTrial(). Tij is the probability that this move would generate state j
	 * from state i, and Tji is the probability that a subsequent call to
	 * doTrial would return to state i from state j.
	 */
	public abstract double lnTrialRatio();

	/**
	 * Returns the log of the limiting-distribution probabilities of states,
	 * ln(Pj/Pi), for the states encountered before (i) and after (j) the most
	 * recent call to doTrial.
	 */
	public abstract double lnProbabilityRatio();

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
	 * move the last time thisTrial was invoked (regardless of whether the move
	 * was accepted). This information usually is not needed, but it is
	 * available in cases where required by objects that perform some activity
	 * while the move is in progress, or need to update after the move is
	 * completed. Such objects can receive notification of the move's completion
	 * by registering with the IntegratorMC as MCMoveEventListeners.
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
	public void setPhase(Phase[] p) {
        if(p.length != phases.length) throw new IllegalArgumentException("Invalid number of phases for MCMove");
		System.arraycopy(p, 0, phases, 0, p.length);
	}

    /**
     * @return the phase(s) on which this move acts.
     */
	public Phase[] getPhase() {
		return phases;
	}

	protected void adjustStepSize() {
		if (nTrials == 0) {
			return;
		}
		nTrialsSum += nTrials;
		nAcceptSum += nAccept;
		if (nTrialsSum != 0)
			acceptanceRatio = (double) nAcceptSum / (double) nTrialsSum;
		if (!tunable) {
			nTrials = nAccept = 0;
			return;
		}
		if (nAccept > (int) (acceptanceTarget * nTrials)) {
			stepSize *= 1.05;
		} else {
			stepSize *= 0.95;
		}
		stepSize = Math.min(stepSize, stepSizeMax);
		stepSize = Math.max(stepSize, stepSizeMin);
		nTrials = 0;
		nAccept = 0;
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
	public final int nominalFrequency() {
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
	 * Zeros all accumulators for acceptance-rate statistics.
	 */
	public void reset() {
		setAcceptanceTarget(acceptanceTarget);
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
	 * Sets the desired rate of acceptance (as a fraction between 0 and 1
	 * inclusive) of trials of this move. IllegalArgumentException is thrown if
	 * given value is outside of acceptable range.
	 */
	public final void setAcceptanceTarget(double target) {
		if (target < 0.0 || target > 1.0) {
			throw new IllegalArgumentException(
					"Acceptance target should be a number between zero and unity");
		}
		nTrialsSum = 0;
		nAcceptSum = 0;
		acceptanceTarget = target;
		nTrials = 0;
		nAccept = 0;
		acceptanceRatio = Double.NaN;
	}

	/**
	 * @return current value of the targeted rate of acceptance for this move.
	 */
	public final double getAcceptanceTarget() {
		return acceptanceTarget;
	}

	public void setStepSize(double step) {
		stepSize = step;
	}

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
	 * @return Returns the temperature.
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * @param temperature
	 *            The temperature to set.
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	public final void setAdjustInterval(int i) {
		adjustInterval = i;
	}

	public final int getAdjustInterval() {
		return adjustInterval;
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

	private double acceptanceRatio, acceptanceTarget;

	protected double stepSize, stepSizeMax, stepSizeMin;

	private int nTrials, nAccept, nTrialsSum, nAcceptSum, adjustInterval;

	protected boolean tunable = true;

	protected final Phase[] phases;

	private String name;

	protected double temperature;

	protected final PotentialMaster potential;

}