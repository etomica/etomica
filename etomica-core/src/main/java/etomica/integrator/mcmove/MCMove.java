/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.atom.iterator.AtomIterator;

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

public abstract class MCMove {

    /**
     * @param potentialMaster the potential master that move can use to calculate energy
     */
	public MCMove(IPotentialMaster potentialMaster) {
        this(potentialMaster, new MCMoveTracker());
    }

    /**
     * @param potentialMaster the potential master that move can use to calculate energy
     */
    public MCMove(IPotentialMaster potentialMaster, MCMoveTracker acceptanceTracker) {
        potential = potentialMaster;
        moveTracker = acceptanceTracker;
        nominalFrequency = 100;
    }

    public void setTracker(MCMoveTracker newTracker) {
        moveTracker = newTracker;
    }

	/**
	 * Method to perform trial move. Returns false if the trial could not be
	 * attempted, for example if there were no molecules in the box and the
	 * trial is designed to displace an atom; returns true otherwise.
	 */
	public abstract boolean doTrial();

	/**
     * Returns the temperature-independent part of the acceptance probability.  
     * The actual acceptance probility is calculated as max(1,A*exp(-B/T))
     * where T is the IntegratorBox's temperature.
	 */
	public abstract double getA();

	/**
	 * Returns the temperature-dependent part of the acceptance probability.  
     * The actual acceptance probility is calculated as max(1,A*exp(-B/T))
     * where T is the IntegratorBox's temperature.
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
     * move caused an atom to be removed from the box, it will not be returned
     * because that Atom is not considered an affected Atom.                    
	 */
	public abstract AtomIterator affectedAtoms(IBox box);

	/**
	 * Returns the change in the energy of the given box that results from the
	 * trial move. Should be called only after lnProbabilityRatio(); returns
	 * Double.NaN if invoked between calls to doTrial and lnProbabilityRatio.
	 */
	public abstract double energyChange(IBox box);

	public MCMoveTracker getTracker() {
        return moveTracker;
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
	
	public String toString() {
	    return getClass().toString();
	}

    protected final IPotentialMaster potential;

    /**
	 * Value giving nominal frequency for performing this move. Default is 100,
	 * but may be given a different value by subclasses.
	 */
	protected int nominalFrequency;
    
    protected MCMoveTracker moveTracker;
}
