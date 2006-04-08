package etomica.integrator.mcmove;

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
        this(potentialMaster, new MCMoveTracker(), nPhases);
    }

    /**
     * @param potentialMaster the potential master that move can use to calculate energy
     * @param nPhases the number of phases on which the move acts.  This is used at
     * construction to size the (final) phases array and cannot be changed.
     */
    public MCMove(PotentialMaster potentialMaster, MCMoveTracker acceptanceTracker, int nPhases) {
        potential = potentialMaster;
        moveTracker = acceptanceTracker;
        nominalFrequency = 100;
        perParticleFrequency = false;
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

	public MCMoveTracker getTracker() {
        return moveTracker;
    }
    
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

    protected Phase phase;
    private String name;
    protected final PotentialMaster potential;

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
    
    protected final MCMoveTracker moveTracker;
}