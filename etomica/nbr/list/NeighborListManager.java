package etomica.nbr.list;

import etomica.action.AtomActionAdapter;
import etomica.action.PhaseImposePbc;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.SpeciesRoot;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorIntervalListener;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.integrator.IntegratorPhase;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.cell.ApiAACell;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentSourceAtomManager;
import etomica.potential.Potential;
import etomica.potential.PotentialArray;
import etomica.util.Arrays;
import etomica.util.Debug;

/**
 * Initiates the process of updating the neighbor lists. Instance is constructed
 * by PotentialMasterNbr constructor. Acts as a listener of the integrator(s),
 * and performs the update at regular intervals upon receiving interval events.
 * Each event causes the manager to loop through all phases acted upon by the
 * integrator (as given by the integrator's getPhase method), and check each
 * atom against any neighbor criteria that apply to it, seeing if it has changed
 * (e.g., moved) in a way that requires its neighbor lists to be updated. When
 * this is found for any atom, all atom neighbor lists are updated via a call to
 * the calculate method of PotentialMasterNbr, passing a
 * PotentialCalculationCellAssign instance as the PotentialCalculation.
 */
public class NeighborListManager implements IntegratorNonintervalListener,
        IntegratorIntervalListener, AgentSource, java.io.Serializable {

    /**
     * Configures instance for use by the given PotentialMaster.
     */
    public NeighborListManager(PotentialMasterList potentialMasterList, double range, 
            PhaseAgentManager agentManager) {
        setUpdateInterval(1);
        iieCount = updateInterval;
        iterator = new AtomIteratorTree();
        iterator.setDoAllNodes(true);
        neighborCheck = new NeighborCheck(this);
        neighborReset = new NeighborReset(this);
        setPriority(200);
        pbcEnforcer = new PhaseImposePbc();
        pbcEnforcer.setApplyToMolecules(false);
        potentialMaster = potentialMasterList;
        cellNbrIterator = new ApiAACell(potentialMaster.getSpace().D(), range, agentManager);
        phaseAgentManager = new PhaseAgentManager(new PhaseAgentSourceAtomManager(this),null);
        phaseAgentManager1Body = new PhaseAgentManager(new PhaseAgentSourceAtomManager(new AtomPotential1ListSource()),null);
        phaseClean = new boolean[0];
    }

    /**
     * Reacts to an integrator INITIALIZE event, preparing the
     * neighbor-list facility. Performs the following actions:
     * <ul>
     * <li>applies periodic boundary conditions, to move all atoms to the
     * central image
     * <li>identifies to each AtomType instance the potentials that apply to
     * atoms of that type
     * <li>assigns each interacting atom to a cell
     * </ul>
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (evt.type() == IntegratorNonintervalEvent.INITIALIZE) {
            Phase phase = ((IntegratorPhase)evt.getSource()).getPhase();
            phaseAgentManager.setRoot((SpeciesRoot)phase.getSpeciesMaster().parentGroup());
            phaseAgentManager1Body.setRoot((SpeciesRoot)phase.getSpeciesMaster().parentGroup());
            updateLists(phase);
            reset(phase);
        }
    }
    
    protected void updateLists(Phase phase) {
        if (phaseClean.length > phase.getIndex() && phaseClean[phase.getIndex()]) {
            // phase is not dirty
            return;
        }
        
        agentManager = (AtomAgentManager)phaseAgentManager.getAgent(phase);
        potentialListManager = (AtomAgentManager)phaseAgentManager1Body.getAgent(phase);

        iterator.setRootAtom(phase.getSpeciesMaster());
        iterator.reset();
        while (iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            int numPotentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials().length;
            ((AtomNeighborLists)agentManager.getAgent(atom)).setCapacity(numPotentials);

            numPotentials = potentialMaster.getIntraPotentials(atom.getType()).getPotentials().length;
            ((AtomPotentialList)potentialListManager.getAgent(atom)).setCapacity(numPotentials);
        }
        
        if (phaseClean.length < phase.getIndex()+1) {
            phaseClean = Arrays.resizeArray(phaseClean,phase.getIndex()+1);
        }
        //mark phase as not dirty
        phaseClean[phase.getIndex()] = true;
    }

    /**
     * Reacts to an interval event, with every updateInterval calls causing
     * this to invoke updateNbrsIfNeeded.
     */
    public void intervalAction(IntegratorIntervalEvent evt) {
        if (--iieCount == 0) {
            updateNbrsIfNeeded(((IntegratorPhase)evt.getSource()));
            iieCount = updateInterval;
        }
    }

    /**
     * For each phase in the array, applies central image, 
     * resets neighbors of all atoms, and sets up all neighbor
     * lists.
     */
    public void reset(Phase phase) {
        for (int j = 0; j < criteriaArray.length; j++) {
            criteriaArray[j].setPhase(phase);
        }
        pbcEnforcer.setPhase(phase);
        pbcEnforcer.actionPerformed();
        neighborSetup(phase);
    }

    /**
     * Checks whether any atom needs neighbor list updating, and if
     * one is found, performs neighbor list updates of all atom 
     * neighbor lists.  Performs this action on all phases acted on
     * by given integrator.
     */
    public void updateNbrsIfNeeded(IntegratorPhase integrator) {
        Phase phase = integrator.getPhase();
        neighborCheck.reset();
        for (int j = 0; j < criteriaArray.length; j++) {
            criteriaArray[j].setPhase(phase);
        }
        iterator.setRootAtom(phase.getSpeciesMaster());
        iterator.allAtoms(neighborCheck);
        if (neighborCheck.needUpdate) {
            if (Debug.ON && Debug.DEBUG_NOW) {
                System.out.println("Updating neighbors");
            }
            if (neighborCheck.unsafe && !quiet) {
                System.err
                        .println("Atoms exceeded the safe neighbor limit");
            }
            pbcEnforcer.setPhase(phase);
            pbcEnforcer.actionPerformed();
            neighborSetup(phase);
            integrator.neighborsUpdated();
        }

    }

    /**
     * Returns the interval for which neighbor update checks are performed.  After receiving
     * this number of interval events, updateNbrsIfNeeded is invoked.
     */
    public int getUpdateInterval() {
        return updateInterval;
    }

    /**
     * Sets the interval for which neighbor update checks are performed.  After receiving
     * this number of interval events, updateNbrsIfNeeded is invoked.
     */
    public void setUpdateInterval(int updateInterval) {
        this.updateInterval = updateInterval;
    }

    /**
     * Adds the given criterion to the list of those held by this manager. This
     * is done so the manager can inform all criteria of the phase in which they
     * are being applied. Does not add the criterion if it was already added to
     * the list.
     */
    public void addCriterion(NeighborCriterion criterion) {
        //mark all phases as unseen (dirty)
        phaseClean = new boolean[0];
        
        boolean found = false;
        for (int i=0; i<criteriaArray.length; i++) {
            if (criteriaArray[i] == criterion) {
                found = true;
                break;
            }
        }
        if (!found) {
            criteriaArray = (NeighborCriterion[]) Arrays.addObject(criteriaArray, criterion);
        }
    }
    
    public NeighborCriterion[] getCriterion(AtomType atomType) {
        return potentialMaster.getRangedPotentials(atomType).getCriteria();
    }

    /**
     * Removes the given criterion from the list of those held by this manager.
     * 
     * @param criterion
     *            Criterion to be removed from the list
     * @return false if the criterion was not previously added to this manager.
     */
    public boolean removeCriterion(NeighborCriterion criterion) {
        int oldLength = criteriaArray.length;
        criteriaArray = (NeighborCriterion[]) Arrays.removeObject(criteriaArray,
                criterion);
        if (oldLength != criteriaArray.length) {
            //we actually cared about this criterion.  mark all phases as unseen (dirty)
            phaseClean = new boolean[0];
            return true;
        }

        return false;
    }

    /**
     * @return Returns the interval-listener priority.
     */
    public int getPriority() {
        return priority;
    }

    /**
     * Sets the interval-listener priority. Default value is 300, which puts
     * this after central-image enforcement.
     * 
     * @param priority
     *            The priority to set.
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
     * @param pbcEnforcer
     *            The pbcEnforcer to set.
     */
    public void setPbcEnforcer(PhaseImposePbc pbcEnforcer) {
        this.pbcEnforcer = pbcEnforcer;
    }

    /**
     * Reassigns all interacting atoms to cells, then loops over all cell-list
     * neighbor pairs, determines for each pair whether a potential applies to it,
     * and if so, puts each in the other's neighbor list.
     * Called by updateNbrsIfNeeded, and by reset.
     * @param phase phase in which neighbor setup is performed.
     */
    protected void neighborSetup(Phase phase) {
        agentManager = (AtomAgentManager)phaseAgentManager.getAgent(phase);
        agentManager1Body = (AtomAgentManager)phaseAgentManager1Body.getAgent(phase);

        iterator.setRootAtom(phase.getSpeciesMaster());
        neighborReset.setNeighborLists(agentManager,agentManager1Body);
        iterator.allAtoms(neighborReset);
        
        potentialMaster.getNbrCellManager(phase).assignCellAll();

        cellNbrIterator.setPhase(phase);
        cellNbrIterator.reset();
        //TODO change looping scheme so getPotentials isn't called for every pair
        //consider doing this by introducing ApiNested interface, with hasNextInner and hasNextOuter methods
        while (cellNbrIterator.hasNext()) {
            AtomPair pair = cellNbrIterator.nextPair();
            Atom atom0 = pair.atom0;
            PotentialArray potentialArray = potentialMaster.getRangedPotentials(atom0.getType());
            Potential[] potentials = potentialArray.getPotentials();
            NeighborCriterion[] criteria = potentialArray.getCriteria();
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() < 2) {
                    continue;
                }
                if (criteria[i].accept(pair)) {
                    ((AtomNeighborLists)agentManager.getAgent(pair.atom0)).addUpNbr(pair.atom1,i);
                    ((AtomNeighborLists)agentManager.getAgent(pair.atom1)).addDownNbr(pair.atom0,
                            potentialMaster.getRangedPotentials(pair.atom1.getType()).getPotentialIndex(potentials[i]));
                }
            }
        }
    }

    /**
     * Sets the interaction range, which affects the cell-list neighbor iteration
     * used to generate candidate neighbors for neighbor listing.
     */
    public void setRange(double d) {
        cellNbrIterator.getNbrCellIterator().setNeighborDistance(d);
    }
    
    public double getRange() {
        return cellNbrIterator.getNbrCellIterator().getNeighborDistance();
    }

    /**
     * @return quiet flag, indicating if unsafe-neighbor conditions should
     *         generate an error message (would not want this if atoms were
     *         inserted in a MC move, for example).
     */
    public boolean isQuiet() {
        return quiet;
    }

    /**
     * Sets the quiet flag, indicating if unsafe-neighbor conditions should generate
     * an error message (would not want this if atoms were inserted in a MC
     * move, for example).
     * 
     * @param quiet
     *            if true, no error will be generated; default is false
     */
    public void setQuiet(boolean quiet) {
        this.quiet = quiet;
    }
    
    public void setPhase(Phase phase) {
        // if agentManagers is null here, it's probably because you haven't added this instance
        // as an IntegratorListener
        agentManager = (AtomAgentManager)phaseAgentManager.getAgent(phase);
        potentialListManager = (AtomAgentManager)phaseAgentManager1Body.getAgent(phase);
    }
    
    public AtomArrayList[] getUpList(Atom atom) {
        return ((AtomNeighborLists)agentManager.getAgent(atom)).getUpList();
    }

    public AtomArrayList[] getDownList(Atom atom) {
        return ((AtomNeighborLists)agentManager.getAgent(atom)).getDownList();
    }

    public AtomPotentialList getPotential1BodyList(Atom atom) {
        return (AtomPotentialList)potentialListManager.getAgent(atom);
    }

    private static final long serialVersionUID = 1L;
    private NeighborCriterion[] criteriaArray = new NeighborCriterion[0];
    private int updateInterval;
    private int iieCount;
    private final AtomIteratorTree iterator;
    private final NeighborCheck neighborCheck;
    private final NeighborReset neighborReset;
    private final ApiAACell cellNbrIterator;
    protected final PotentialMasterList potentialMaster;
    private int priority;
    private PhaseImposePbc pbcEnforcer;
    private boolean quiet;
    private AtomAgentManager potentialListManager;
    private AtomAgentManager agentManager;
    private AtomAgentManager agentManager1Body;
    private final PhaseAgentManager phaseAgentManager;
    private final PhaseAgentManager phaseAgentManager1Body;
    private boolean[] phaseClean;

    /**
     * Atom action class that checks if any criteria indicate that the given
     * atom needs to update its neighbor list.
     */
    private static class NeighborCheck extends AtomActionAdapter {

        private static final long serialVersionUID = 1L;
        protected boolean needUpdate = false, unsafe = false;
        private NeighborListManager neighborListManager;

        public NeighborCheck(NeighborListManager manager) {
            neighborListManager = manager;
        }
        
        public void actionPerformed(Atom atom) {
            final NeighborCriterion[] criterion = neighborListManager.getCriterion(atom.getType());
            for (int i = 0; i < criterion.length; i++) {
                if (criterion[i].needUpdate(atom)) {
                    needUpdate = true;
                    if (criterion[i].unsafe()) {
                        if (Debug.ON && Debug.DEBUG_NOW) {
                            System.out.println("atom " + atom
                                    + " exceeded safe limit");
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
     * Atom action class that clears neighbor list of given atom and loops
     * through all neighbor criteria applying to atom (as given by its type),
     * and resets the criteria as it applies to the atom (e.g., sets its
     * previous-position vector to its current position).
     */
    private static class NeighborReset extends AtomActionAdapter {
        private static final long serialVersionUID = 1L;

        public NeighborReset(NeighborListManager manager) {
            neighborListManager = manager;
        }
        
        public void actionPerformed(Atom atom) {
            //TODO consider removing this check, for perf improvement
            if (atom.getType().getDepth() < 3) {
                return;//don't want SpeciesMaster or SpeciesAgents
            }
            final NeighborCriterion[] criterion = neighborListManager.getCriterion(atom.getType());
            ((AtomNeighborLists)agentManager.getAgent(atom)).clearNbrs();
            for (int i = 0; i < criterion.length; i++) {
                criterion[i].reset(atom);
            }

            PotentialArray potentialArray = neighborListManager.potentialMaster.getRangedPotentials(atom.getType());
            Potential[] potentials = potentialArray.getPotentials();
            NeighborCriterion[] criteria = potentialArray.getCriteria();

            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() != 1) {
                    continue;
                }
                ((AtomPotentialList)agentManager1Body.getAgent(atom)).setIsInteracting(criteria[i].accept(atom),i);
            }
        }
        
        public void setNeighborLists(AtomAgentManager newAgentManager, AtomAgentManager newAgentManager1Body) {
            agentManager = newAgentManager;
            agentManager1Body = newAgentManager1Body;
        }
        
        private NeighborListManager neighborListManager;
        private AtomAgentManager agentManager;
        private AtomAgentManager agentManager1Body;
    }
    
    public Class getAgentClass() {
        return AtomNeighborLists.class;
    }
    
    public Object makeAgent(Atom atom) {
        AtomNeighborLists lists = new AtomNeighborLists();
        int numPotentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials().length;
        lists.setCapacity(numPotentials);
        return lists;
    }
    
    public void releaseAgent(Object agent, Atom atom) {
        ((AtomNeighborLists)agent).clearNbrs();
    }
    
    public class AtomPotential1ListSource implements AtomAgentManager.AgentSource, java.io.Serializable {
        private static final long serialVersionUID = 1L;
        public Class getAgentClass() {
            return AtomPotentialList.class;
        }
        public void releaseAgent(Object obj, Atom atom) {}
        public Object makeAgent(Atom atom) {
            AtomPotentialList lists = new AtomPotentialList();
            int numPotentials = potentialMaster.getIntraPotentials(atom.getType()).getPotentials().length;
            lists.setCapacity(numPotentials);
            return lists;
        }
    }
}
