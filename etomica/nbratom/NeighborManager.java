/*
 * History
 * Created on Sep 22, 2004 by kofke
 */
package etomica.nbratom;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomSet;
import etomica.Debug;
import etomica.Integrator;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Potential;
import etomica.Integrator.IntervalEvent;
import etomica.Integrator.IntervalListener;
import etomica.action.AtomsetActionAdapter;
import etomica.action.PhaseImposePbc;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.nbr.NeighborCriterion;
import etomica.nbratom.cell.ApiAACell;
import etomica.nbratom.cell.AtomIteratorCell;
import etomica.nbratom.cell.NeighborCellManager;
import etomica.nbratom.cell.PotentialCalculationCellAssign;
import etomica.potential.Potential2;
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
public class NeighborManager implements IntervalListener {

	/**
	 * Configures instance for use by the given PotentialMaster.
	 */    
	public NeighborManager(PotentialMasterNbr potentialMaster) {
		super();
		setUpdateInterval(1);
		iieCount = updateInterval;
		iterator = new AtomIteratorTree();
		iterator.setDoAllNodes(true);
		neighborCheck = new NeighborCheck();
		neighborReset = new NeighborReset();
        setPriority(200);
        pbcEnforcer = new PhaseImposePbc();
        pbcEnforcer.setApplyToMolecules(false);
        this.potentialMaster = potentialMaster;
        cellNbrIterator = new ApiAACell(potentialMaster.getSpace().D());
        atomIterator = new AtomIteratorCell(potentialMaster.getSpace().D());
	}

	/* (non-Javadoc)
	 * @see etomica.Integrator.IntervalListener#intervalAction(etomica.Integrator.IntervalEvent)
	 */
	public void intervalAction(IntervalEvent evt) {
		if(evt.type() == IntervalEvent.START || evt.type() == IntervalEvent.INITIALIZE) {
            Phase[] phases = ((Integrator)evt.getSource()).getPhase();
            IteratorDirective idUp = new IteratorDirective();
            for (int i=0; i<phases.length; i++) {
                pbcEnforcer.setPhase(phases[i]);
                pbcEnforcer.actionPerformed();
                PotentialCalculationCellAssign pc = new PotentialCalculationCellAssign(potentialMaster.getNbrCellManager(phases[i]));
                potentialMaster.calculate(phases[i],idUp,pc);
            }
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
            pbcEnforcer.setPhase(phase[i]);
            pbcEnforcer.actionPerformed();
			iterator.setRoot(phase[i].speciesMaster());
			iterator.allAtoms(neighborReset);
			neighborSetup(phase[i]);
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
                neighborSetup(phase[i]);
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

    public void neighborSetup(Phase phase) {
        NeighborCellManager nbrManager = potentialMaster.getNbrCellManager(phase);
        atomIterator.setPhase(phase);
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            nbrManager.assignCell(atomIterator.nextAtom());
        }
        
        cellNbrIterator.setPhase(phase);
        cellNbrIterator.reset();
        //TODO change looping scheme so getPotentials isn't called for every pair
        while (cellNbrIterator.hasNext()) {
            AtomPair pair = cellNbrIterator.nextPair();
            Atom atom0 = pair.atom0;
            Potential[] potentials = atom0.type.getNbrManagerAgent().getPotentials();
            for (int i=0; i<potentials.length; i++) {
                if (((Potential2)potentials[i]).getCriterion().accept(pair)) {
                    ((AtomSequencerNbr)pair.atom0.seq).addUpNbr(pair.atom1, potentials[i]);
                    ((AtomSequencerNbr)pair.atom1.seq).addDownNbr(pair.atom0, potentials[i]);
                }
            }
        }
    }
    
    public void setRange(double d) {
        cellNbrIterator.getNbrCellIterator().setRange(d);
    }
    
	private NeighborCriterion[] criteria = new NeighborCriterion[0];
	private int updateInterval;
	private int iieCount;
	private final AtomIteratorTree iterator;
	private final NeighborCheck neighborCheck;
	private final NeighborReset neighborReset;
    private final ApiAACell cellNbrIterator;
    private final AtomIteratorCell atomIterator;
    private final PotentialMasterNbr potentialMaster;
	private int priority;
    private PhaseImposePbc pbcEnforcer;
    
	private static class NeighborCheck extends AtomsetActionAdapter {
		private boolean needUpdate = false, unsafe = false;
		public void actionPerformed(AtomSet atom) {
			final NeighborCriterion[] criterion = ((Atom)atom).type.getNbrManagerAgent().getCriterion();
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
		
		public void reset() {
			needUpdate = false;
			unsafe = false;
		}

	}
	
	private class NeighborReset extends AtomsetActionAdapter {
		public void actionPerformed(AtomSet atom) {
            //FIXME changes to depth might make this wrong
            if(((Atom)atom).type.getDepth() < 2) return;//don't want SpeciesMaster or SpeciesAgents
			final NeighborCriterion[] criterion = ((Atom)atom).type.getNbrManagerAgent().getCriterion();
			((AtomSequencerNbr)((Atom)atom).seq).clearNbrs();
			for (int i=0; i<criterion.length; i++) {
				criterion[i].reset((Atom)atom);
			}
		}
	}
}
