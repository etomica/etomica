/*
 * History
 * Created on Sep 22, 2004 by kofke
 */
package etomica.nbr;

import etomica.Atom;
import etomica.AtomIteratorTree;
import etomica.Debug;
import etomica.Integrator;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Space;
import etomica.Integrator.IntervalEvent;
import etomica.Integrator.IntervalListener;
import etomica.action.AtomsetActionAdapter;
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
 * PotentialCalculationNbrSetup instance as the PotentialCalculation.  
 */
public class NeighborManager implements IntervalListener {

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
	}

	/* (non-Javadoc)
	 * @see etomica.Integrator.IntervalListener#intervalAction(etomica.Integrator.IntervalEvent)
	 */
	public void intervalAction(IntervalEvent evt) {
		if(evt.type() == IntervalEvent.START) {
			reset(((Integrator)evt.getSource()).getPhase());
		} else if(evt.type() == IntervalEvent.INTERVAL) {
			if (--iieCount == 0) {
				updateNbrsIfNeeded((Integrator)evt.getSource());
				iieCount = updateInterval;
			}
		}
	}

	/**
	 * @param 
	 */
	public void reset(Phase[] phase) {
		for (int i=0; i<phase.length; i++) {
			for (int j=0; j<criteria.length; j++) {
				criteria[j].setPhase(phase[i]);
			}
			boundary = phase[i].boundary();
			iterator.setRoot(phase[i].speciesMaster());
			iterator.allAtoms(neighborReset);
			potentialMaster.calculate(phase[i],id,potentialCalculationNbrSetup);
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
				boundary = phase[i].boundary();
				iterator.allAtoms(neighborReset);
				potentialMaster.calculate(phase[i],id,potentialCalculationNbrSetup);
				resetIntegrator = true;
			}
		}
        //TODO consider a reset(Phase) method for integrator to reset relative to just the affected phase
        if(resetIntegrator) integrator.neighborsUpdated();
	}
	
	public int getUpdateInterval() {
		return updateInterval;
	}
	
	public void setUpdateInterval(int updateInterval) {
		this.updateInterval = updateInterval;
	}
	
    /**
     * Adds the given criterion to the list of those held by this manager.
     * Does not check whether criterion already was added to list.
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


	private NeighborCriterion[] criteria = new NeighborCriterion[0];
	private int updateInterval;
	private int iieCount;
	private Space.Boundary boundary;
	private final PotentialMasterNbr potentialMaster;
	private final AtomIteratorTree iterator;
	private final NeighborCheck neighborCheck;
	private final NeighborReset neighborReset;
	private final IteratorDirective id = new IteratorDirective();
	private final PotentialCalculationNbrSetup potentialCalculationNbrSetup = new PotentialCalculationNbrSetup();
	private int priority;
    
	private static class NeighborCheck extends AtomsetActionAdapter {
		private boolean needUpdate = false, unsafe = false;
		public void actionPerformed(Atom[] atom) {
			NeighborCriterion criterion = atom[0].type.getNbrManagerAgent().getCriterion();
			if (criterion != null && criterion.needUpdate(atom[0])) {
				needUpdate = true;
				if (criterion.unsafe()) {
					if (Debug.DEBUG_NOW) {
						System.out.println("atom "+atom[0]+" exceeded safe limit");
					}
					unsafe = true;
				}
			}
		}
		
		public void reset() {
			needUpdate = false;
			unsafe = false;
		}

	}
	
	private class NeighborReset extends AtomsetActionAdapter {
		public void actionPerformed(Atom[] atom) {
            if(atom[0].node.depth() < 2) return;//don't want SpeciesMaster or SpeciesAgents
			NeighborCriterion criterion = atom[0].type.getNbrManagerAgent().getCriterion();
			((AtomSequencerNbr)atom[0].seq).clearNbrs();
			if (criterion != null) {
				boundary.centralImage(atom[0].coord);
				criterion.reset(atom[0]);
			}
		}
	}
}
