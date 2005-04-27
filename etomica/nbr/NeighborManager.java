/*
 * History
 * Created on Sep 22, 2004 by kofke
 */
package etomica.nbr;

import etomica.Atom;
import etomica.AtomSet;
import etomica.Debug;
import etomica.Integrator;
import etomica.IntegratorEvent;
import etomica.IntegratorIntervalEvent;
import etomica.IntegratorIntervalListener;
import etomica.IntegratorListener;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.action.AtomsetActionAdapter;
import etomica.action.PhaseImposePbc;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.utility.Arrays;

/**
 * Initiates the process of updating the neighbor lists.   
 * Instance is constructed by PotentialMasterNbr constructor.
 * Acts as a listener of the integrator(s), and performs the update at regular intervals
 * upon receiving interval events.  Each event causes the manager to loop through
 * all phases acted upon by the integrator (as given by the integrator's getPhase
 * method), and check each atom against any neighbor criteria that apply to it, seeing
 * if it has changed (e.g., moved) in a way that requires its neighbor lists to be
 * updated.  When this is found for any atom, all atom neighbor lists are updated
 * via a call to the calculate method of PotentialMasterNbr, passing a 
 * PotentialCalculationCellAssign instance as the PotentialCalculation.  
 */
public class NeighborManager implements IntegratorListener, IntegratorIntervalListener {

	/**
	 * Configures instance for use by the given PotentialMaster.
	 */    
	public NeighborManager(PotentialMasterNbr potentialMaster) {
		super();
		this.potentialMaster = potentialMaster;
		setUpdateInterval(1);
		iieCount = updateInterval;
		iterator = new AtomIteratorTree();
		iterator.setDoAllNodes(true);
		neighborCheck = new NeighborCheck();
		neighborReset = new NeighborReset();
        setPriority(200);
        pbcEnforcer = new PhaseImposePbc();
        pbcEnforcer.setApplyToMolecules(true);
	}

	/* (non-Javadoc)
	 * @see etomica.Integrator.IntervalListener#intervalAction(etomica.Integrator.IntervalEvent)
	 */
    public void integratorAction(IntegratorEvent evt) {
		if((evt.type().mask & (IntegratorEvent.START.mask | IntegratorEvent.INITIALIZE.mask)) != 0) {
			reset(((Integrator)evt.getSource()).getPhase());
        }
    }
    
    public void intervalAction(IntegratorIntervalEvent evt) {
		if (--iieCount == 0) {
			updateNbrsIfNeeded((Integrator)evt.getSource());
			iieCount = updateInterval;
		}
	}

	/**
	 * Method to initialize the neighbor list for the given phases.
     * For each given phase, applies pbc, clears neighbor lists of 
     * all atoms, and performs setup of neighbor lists for current
     * configuration. 
	 */
	public void reset(Phase[] phase) {
		for (int i=0; i<phase.length; i++) {
			for (int j=0; j<criteria.length; j++) {
				criteria[j].setPhase(phase[i]);
			}
            pbcEnforcer.setPhase(phase[i]);
            pbcEnforcer.actionPerformed();
			iterator.setRoot(phase[i].speciesMaster());
			iterator.allAtoms(neighborReset);
			potentialMaster.calculate(phase[i], id, potentialCalculationNbrSetup);
		}
	}

    
	public void updateNbrsIfNeeded(Integrator integrator) {
        boolean resetIntegrator = false;
        Phase[] phase = integrator.getPhase();
		for (int i=0; i<phase.length; i++) {
			neighborCheck.reset();
			for (int j=0; j<criteria.length; j++) {
				criteria[j].setPhase(phase[i]);
			}
			iterator.setRoot(phase[i].speciesMaster());
			iterator.allAtoms(neighborCheck);
			if (neighborCheck.needUpdate) {
				if (Debug.ON && Debug.DEBUG_NOW) {
					System.out.println("Updating neighbors");
				}
				if (neighborCheck.unsafe) {
					System.err.println("Atoms exceeded the safe neighbor limit");
				}
                pbcEnforcer.setPhase(phase[i]);
                pbcEnforcer.actionPerformed();
				iterator.allAtoms(neighborReset);
				potentialMaster.calculate(phase[i],id,potentialCalculationNbrSetup);
				resetIntegrator = true;
			}
		}
        //TODO consider a reset(Phase) method for integrator to reset relative to just the affected phase
        if(resetIntegrator) integrator.neighborsUpdated();
	}
	
    /**
     * Returns the number of interval events received before the interval action
     * (i.e., check for need to update neighbor lists) is performed.  Default is 1.
     */
	public int getUpdateInterval() {
		return updateInterval;
	}
	
    /**
     * Sets the number of interval events received before the interval action
     * (i.e., check for need to update neighbor lists) is performed.
     */
	public void setUpdateInterval(int updateInterval) {
		this.updateInterval = updateInterval;
	}
	
    /**
     * Adds the given criterion to the list of those held by this manager.
     * This is done so the manager can inform all criteria of the phase in which
     * they are being applied. Does not check whether criterion already was added to list.
     */
	public void addCriterion(NeighborCriterion criterion) {
        criteria = (NeighborCriterion[])Arrays.addObject(criteria, criterion);
	}

    /**
     * Removes the given criterion from the list of those held by this manager.
     * @param criterion Criterion to be removed from the list
     * @return false if the criterion was not previously added to this manager.
     */
    public boolean removeCriterion(NeighborCriterion criterion) {
        int oldLength = criteria.length;
        criteria = (NeighborCriterion[])Arrays.removeObject(criteria, criterion);
        return (oldLength != criteria.length);
    }
    
    /**
     * @return Returns the interval-listener priority.
     */
    public int getPriority() {
        return priority;
    }
    /**
     * Sets the interval-listener priority.  Default value is 300, which
     * puts this after central-image enforcement.
     * @param priority The priority to set.
     */
    public void setPriority(int priority) {
        this.priority = priority;
    }

    /**
     * @return Returns the pbcEnforcer.
     */
    public PhaseImposePbc getPbcEnforcer() {
        return pbcEnforcer;
    }
    /**
     * @param pbcEnforcer The pbcEnforcer to set.
     */
    public void setPbcEnforcer(PhaseImposePbc pbcEnforcer) {
        this.pbcEnforcer = pbcEnforcer;
    }

	private NeighborCriterion[] criteria = new NeighborCriterion[0];
	private int updateInterval;
	private int iieCount;
	private final PotentialMasterNbr potentialMaster;
	private final AtomIteratorTree iterator;
	private final NeighborCheck neighborCheck;
	private final NeighborReset neighborReset;
	private final IteratorDirective id = new IteratorDirective();
	private final PotentialCalculationNbrSetup potentialCalculationNbrSetup = new PotentialCalculationNbrSetup();
	private int priority;
    private PhaseImposePbc pbcEnforcer;
    
    /**
     * Atom action class that checks if any criteria indicate that the given atom
     * needs to update its neighbor list.
     */
	private static class NeighborCheck extends AtomsetActionAdapter {
		private boolean needUpdate = false, unsafe = false;
		public void actionPerformed(AtomSet atom) {
			NeighborCriterion[] criterion = ((Atom)atom).type.getNbrManagerAgent().getCriterion();
            for (int i=0; i<criterion.length; i++) {
                if (criterion[i].needUpdate((Atom)atom)) {
                    needUpdate = true;
                    if (criterion[i].unsafe()) {
                        if (Debug.DEBUG_NOW) {
                            System.out.println("atom "+atom+" exceeded safe limit");
                        }
                        unsafe = true;
                    }
                }
            }
		}
		
        /**
         * Sets class to condition that indicates that no atoms need to have
         * their neighbor list updated.
         */
		public void reset() {
			needUpdate = false;
			unsafe = false;
		}

	}
	
    /**
     * Atom action class that clears neighbor list of given atom and
     * loops through all neighbor criteria applying to atom (as given 
     * by its type), and resets the criteria as it applies to the atom 
     * (e.g., sets its previous-position vector to its current position).
     */
	private class NeighborReset extends AtomsetActionAdapter {
		public void actionPerformed(AtomSet atom) {
            if(((Atom)atom).type.getDepth() < 3) return;//don't want SpeciesMaster or SpeciesAgents
			final NeighborCriterion[] criterion = ((Atom)atom).type.getNbrManagerAgent().getCriterion();
			((AtomSequencerNbr)((Atom)atom).seq).clearNbrs();
			for (int i=0; i<criterion.length; i++) {
				criterion[i].reset((Atom)atom);
			}
		}
	}
}
