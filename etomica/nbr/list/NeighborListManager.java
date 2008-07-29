package etomica.nbr.list;

import java.io.Serializable;

import etomica.action.AtomAction;
import etomica.action.BoxImposePbc;
import etomica.api.IAction;
import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IIntegratorNonintervalListener;
import etomica.api.IPotential;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.ApiAACell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.PotentialArray;
import etomica.space.ISpace;
import etomica.util.Debug;

/**
 * Initiates the process of updating the neighbor lists. Instance is constructed
 * by PotentialMasterNbr constructor. Acts as a listener of the integrator(s),
 * and performs the update at regular intervals upon receiving interval events.
 * Each event causes the manager to loop through all boxs acted upon by the
 * integrator (as given by the integrator's getBox method), and check each
 * atom against any neighbor criteria that apply to it, seeing if it has changed
 * (e.g., moved) in a way that requires its neighbor lists to be updated. When
 * this is found for any atom, all atom neighbor lists are updated via a call to
 * the calculate method of PotentialMasterNbr, passing a
 * PotentialCalculationCellAssign instance as the PotentialCalculation.
 */
public class NeighborListManager implements IIntegratorNonintervalListener,
        IAction, AgentSource, Serializable {

    /**
     * Configures instance for use by the given PotentialMaster.
     */
    public NeighborListManager(PotentialMasterList potentialMasterList, double range, 
            IBox box, ISpace space) {
        setUpdateInterval(1);
        this.box = box;
        iieCount = updateInterval;
        neighborCheck = new NeighborCheck(this);
        setPriority(200);
        pbcEnforcer = new BoxImposePbc(space);
        pbcEnforcer.setBox(box);
        pbcEnforcer.setApplyToMolecules(false);
        potentialMaster = potentialMasterList;
        cellNbrIterator = new ApiAACell(space.D(), range, box);
        cell1ANbrIterator = new Api1ACell(space.D(), range, potentialMasterList.getCellAgentManager());
        agentManager2Body = new AtomLeafAgentManager(this, box);
        agentManager1Body = new AtomLeafAgentManager(new AtomPotential1ListSource(potentialMasterList), box);
        neighborReset = new NeighborReset(this, agentManager2Body, agentManager1Body);
        boxEvent = new BoxEventNeighborsUpdated(box);
        initialized = false;
        doApplyPBC = true;
    }

    public void setDoApplyPBC(boolean newDoApplyPBC) {
        doApplyPBC = newDoApplyPBC;
    }

    public boolean getDoApplyPBC() {
        return doApplyPBC;
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
        if (evt.type() == IntegratorNonintervalEvent.RESET) {
            reset();
        }
    }
    
    public void updateLists() {
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int j=0; j<nLeaf; j++) {
            IAtom atom = leafList.getAtom(j);
            IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials();

            ((AtomNeighborLists)agentManager2Body.getAgent((IAtomLeaf)atom)).setCapacity(potentials.length);
            ((AtomPotentialList)agentManager1Body.getAgent((IAtomLeaf)atom)).setCapacity(potentials.length);
        }
    }

    /**
     * Reacts to an interval event, with every updateInterval calls causing
     * this to invoke updateNbrsIfNeeded.
     */
    public void actionPerformed() {
        if (--iieCount == 0) {
            updateNbrsIfNeeded();
            iieCount = updateInterval;
        }
    }

    /**
     * For each box in the array, applies central image, 
     * resets neighbors of all atoms, and sets up all neighbor
     * lists.
     */
    public void reset() {
        // the NeighborCellManager might not have existed during construction
        // so we couldn't se the lattice.  It better exist now.
        cellNbrIterator.setLattice(potentialMaster.getNbrCellManager(box).getLattice());

        NeighborCriterion[] criteriaArray = potentialMaster.getNeighborCriteria();
        if (oldCriteria != criteriaArray) {
            // if the array of criteria is different, a potential was added or
            // removed and we need to update our lists
            updateLists();
            oldCriteria = criteriaArray;
        }
        for (int j = 0; j < criteriaArray.length; j++) {
            criteriaArray[j].setBox(box);
        }
        if (doApplyPBC) {
            pbcEnforcer.actionPerformed();
        }
        neighborSetup();
        iieCount = updateInterval;
    }

    /**
     * Checks whether any atom needs neighbor list updating, and if
     * one is found, performs neighbor list updates of all atom 
     * neighbor lists.  Performs this action on all boxs acted on
     * by given integrator.
     */
    public void updateNbrsIfNeeded() {
        neighborCheck.reset();
        NeighborCriterion[] criteriaArray = potentialMaster.getNeighborCriteria();
        if (oldCriteria != criteriaArray) {
            // if the array of criteria is different, a potential was added or
            // removed and we need to update our lists
            updateLists();
            oldCriteria = criteriaArray;
        }
        for (int j = 0; j < criteriaArray.length; j++) {
            criteriaArray[j].setBox(box);
        }

        
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int j=0; j<nLeaf; j++) {
            neighborCheck.actionPerformed(leafList.getAtom(j));
        }
        
        if (neighborCheck.needUpdate) {
            if (Debug.ON && Debug.DEBUG_NOW) {
                System.out.println("Updating neighbors");
            }
            if (neighborCheck.unsafe && !quiet) {
                System.err
                        .println("Atoms exceeded the safe neighbor limit");
            }
            if (doApplyPBC) {
                pbcEnforcer.actionPerformed();
            }
            neighborSetup();
            box.getEventManager().fireEvent(boxEvent);
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
        iieCount = updateInterval;
    }

    public NeighborCriterion[] getCriterion(IAtomType atomType) {
        return potentialMaster.getRangedPotentials(atomType).getCriteria();
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
    public BoxImposePbc getPbcEnforcer() {
        return pbcEnforcer;
    }

    /**
     * @param pbcEnforcer
     *            The pbcEnforcer to set.
     */
    public void setPbcEnforcer(BoxImposePbc pbcEnforcer) {
        this.pbcEnforcer = pbcEnforcer;
    }

    /**
     * Reassigns all interacting atoms to cells, then loops over all cell-list
     * neighbor pairs, determines for each pair whether a potential applies to it,
     * and if so, puts each in the other's neighbor list.
     * Called by updateNbrsIfNeeded, and by reset.
     * @param box box in which neighbor setup is performed.
     */
    protected void neighborSetup() {

        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int j=0; j<nLeaf; j++) {
            neighborReset.actionPerformed(leafList.getAtom(j));
        }
        
        NeighborCellManager cellManager = potentialMaster.getNbrCellManager(box);
        cellManager.setDoApplyPBC(!doApplyPBC);
        cellManager.assignCellAll();

        cellNbrIterator.reset();
        //TODO change looping scheme so getPotentials isn't called for every pair
        //consider doing this by introducing ApiNested interface, with hasNextInner and hasNextOuter methods
        for (IAtomSet pair = cellNbrIterator.nextPair(); pair != null;
             pair = cellNbrIterator.nextPair()) {
            IAtom atom0 = pair.getAtom(0);
            IAtom atom1 = pair.getAtom(1);
            PotentialArray potentialArray = potentialMaster.getRangedPotentials(atom0.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            NeighborCriterion[] criteria = potentialArray.getCriteria();
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() < 2) {
                    continue;
                }
                if (criteria[i].accept(pair)) {
                    ((AtomNeighborLists)agentManager2Body.getAgent((IAtomLeaf)atom0)).addUpNbr(atom1,i);
                    ((AtomNeighborLists)agentManager2Body.getAgent((IAtomLeaf)atom1)).addDownNbr(atom0,
                            potentialMaster.getRangedPotentials(atom1.getType()).getPotentialIndex(potentials[i]));
                }
            }
        }
        initialized = true;
    }

    /**
     * Constructs neighbor lists for the given atom
     */
    public void addAtomNotify(IAtomLeaf atom) {
        if (!initialized) {
            // the simulation hasn't started yet.  just wait for neighborSetup
            // to get called.  It can do everything at once and can be sure
            // things are set up properly (like calling setBox on the criteria).
            return;
        }
        if (agentManager2Body.getAgent(atom) == null) {
            // we're getting called before our own makeAgent (we have no
            // control over order here).  make the agent now and then use it
            // again from makeAgent (ugh).  this depends on AtomAgentManager
            // nulling out agents for removed atoms.
            agentManager2Body.setAgent(atom, makeAgent(atom));
        }
        cell1ANbrIterator.setBox(box);
        cell1ANbrIterator.setTarget(atom);
        cell1ANbrIterator.reset();
        for (IAtomSet pair = cell1ANbrIterator.next(); pair != null;
             pair = cell1ANbrIterator.next()) {
            IAtom atom0 = pair.getAtom(0);
            IAtom atom1 = pair.getAtom(1);
            PotentialArray potentialArray = potentialMaster.getRangedPotentials(atom0.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            NeighborCriterion[] criteria = potentialArray.getCriteria();
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() < 2) {
                    continue;
                }
                if (criteria[i].accept(pair)) {
                    ((AtomNeighborLists)agentManager2Body.getAgent((IAtomLeaf)atom0)).addUpNbr(atom1,i);
                    ((AtomNeighborLists)agentManager2Body.getAgent((IAtomLeaf)atom1)).addDownNbr(atom0,
                            potentialMaster.getRangedPotentials(atom1.getType()).getPotentialIndex(potentials[i]));
                }
            }
        }
    }

    /**
     * Sets the interaction range, which affects the cell-list neighbor iteration
     * used to generate candidate neighbors for neighbor listing.
     */
    public void setRange(double d) {
        cell1ANbrIterator.getNbrCellIterator().setNeighborDistance(d);
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
    
    public IAtomSet[] getUpList(IAtom atom) {
        return ((AtomNeighborLists)agentManager2Body.getAgent((IAtomLeaf)atom)).getUpList();
    }

    public IAtomSet[] getDownList(IAtom atom) {
        return ((AtomNeighborLists)agentManager2Body.getAgent((IAtomLeaf)atom)).getDownList();
    }

    public AtomPotentialList getPotential1BodyList(IAtom atom) {
        return (AtomPotentialList)agentManager1Body.getAgent((IAtomLeaf)atom);
    }
    
    public void dispose() {
        agentManager1Body.dispose();
        agentManager2Body.dispose();
    }

    private static final long serialVersionUID = 1L;
    private int updateInterval;
    private int iieCount;
    private final NeighborCheck neighborCheck;
    private final NeighborReset neighborReset;
    private final ApiAACell cellNbrIterator;
    private final Api1ACell cell1ANbrIterator;
    protected final PotentialMasterList potentialMaster;
    private int priority;
    private BoxImposePbc pbcEnforcer;
    private boolean quiet;
    private final AtomLeafAgentManager agentManager2Body;
    private final AtomLeafAgentManager agentManager1Body;
    protected IBox box;
    private NeighborCriterion[] oldCriteria;
    protected final BoxEventNeighborsUpdated boxEvent;
    protected boolean initialized;
    protected boolean doApplyPBC;

    /**
     * Atom action class that checks if any criteria indicate that the given
     * atom needs to update its neighbor list.
     */
    private static class NeighborCheck implements AtomAction, Serializable {

        private static final long serialVersionUID = 1L;
        protected boolean needUpdate = false, unsafe = false;
        private NeighborListManager neighborListManager;

        public NeighborCheck(NeighborListManager manager) {
            neighborListManager = manager;
        }
        
        public void actionPerformed(IAtom atom) {
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
    private static class NeighborReset implements AtomAction, Serializable {
        private static final long serialVersionUID = 1L;

        public NeighborReset(NeighborListManager manager, AtomLeafAgentManager agentManager2Body,
                AtomLeafAgentManager agentManager1Body) {
            neighborListManager = manager;
            this.agentManager2Body = agentManager2Body;
            this.agentManager1Body = agentManager1Body;
            atomSetSinglet = new AtomSetSinglet();
        }
        
        public void actionPerformed(IAtom atom) {
            final NeighborCriterion[] criterion = neighborListManager.getCriterion(atom.getType());
            ((AtomNeighborLists)agentManager2Body.getAgent((IAtomLeaf)atom)).clearNbrs();
            for (int i = 0; i < criterion.length; i++) {
                criterion[i].reset(atom);
            }

            PotentialArray potentialArray = neighborListManager.potentialMaster.getRangedPotentials(atom.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            NeighborCriterion[] criteria = potentialArray.getCriteria();

            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() != 1) {
                    continue;
                }
                atomSetSinglet.atom = atom;
                ((AtomPotentialList)agentManager1Body.getAgent((IAtomLeaf)atom)).setIsInteracting(criteria[i].accept(atomSetSinglet),i);
            }
        }
        
        private NeighborListManager neighborListManager;
        private AtomLeafAgentManager agentManager2Body;
        private AtomLeafAgentManager agentManager1Body;
        protected final AtomSetSinglet atomSetSinglet;
    }
    
    public Class getAgentClass() {
        return AtomNeighborLists.class;
    }
    
    public Object makeAgent(IAtomLeaf atom) {
        if (initialized) {
            Object oldAgent = agentManager2Body.getAgent(atom);
            if (oldAgent != null) {
                // NeighborCellManager got notified first and we already made the
                // agent (and found the neighbors!).  Return that now.
                return oldAgent;
            }
        }
        AtomNeighborLists lists = new AtomNeighborLists();
        IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials();
        lists.setCapacity(potentials.length);
        return lists;
    }
    
    public void releaseAgent(Object agent, IAtomLeaf atom) {
        // we need to remove this atom from the neighbor lists of its neighbors.
        AtomNeighborLists nbrLists = (AtomNeighborLists)agent;
        IAtomSet[] upDnLists = nbrLists.getUpList();
        for (int i=0; i<upDnLists.length; i++) {
            int nNbrs = upDnLists[i].getAtomCount();
            for (int j=0; j<nNbrs; j++) {
                IAtomLeaf jAtom = (IAtomLeaf)upDnLists[i].getAtom(j);
                AtomNeighborLists jNbrLists = (AtomNeighborLists)agentManager2Body.getAgent(jAtom);
                AtomArrayList[] jDnLists = jNbrLists.downList;
                for (int k=0; k<jDnLists.length; k++) {
                    int idx = jDnLists[k].indexOf(atom);
                    if (idx > -1) {
                        jDnLists[k].removeAndReplace(idx);
                    }
                }
            }
        }
        upDnLists = nbrLists.getDownList();
        for (int i=0; i<upDnLists.length; i++) {
            int nNbrs = upDnLists[i].getAtomCount();
            for (int j=0; j<nNbrs; j++) {
                IAtomLeaf jAtom = (IAtomLeaf)upDnLists[i].getAtom(j);
                AtomNeighborLists jNbrLists = (AtomNeighborLists)agentManager2Body.getAgent(jAtom);
                AtomArrayList[] jUpLists = jNbrLists.upList;
                for (int k=0; k<jUpLists.length; k++) {
                    int idx = jUpLists[k].indexOf(atom);
                    if (idx > -1) {
                        jUpLists[k].removeAndReplace(idx);
                    }
                }
            }
        }
        nbrLists.clearNbrs();
    }

    public static class AtomPotential1ListSource implements AtomLeafAgentManager.AgentSource, java.io.Serializable {
        private static final long serialVersionUID = 2L;
        protected final PotentialMasterList potentialMaster;

        public AtomPotential1ListSource(PotentialMasterList potentialMaster) {
            this.potentialMaster = potentialMaster;
        }
        
        public Class getAgentClass() {
            return AtomPotentialList.class;
        }
        public void releaseAgent(Object obj, IAtomLeaf atom) {}
        public Object makeAgent(IAtomLeaf atom) {
            AtomPotentialList lists = new AtomPotentialList();
            IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials();
            
            lists.setCapacity(potentials.length);
            return lists;
        }
    }
}
