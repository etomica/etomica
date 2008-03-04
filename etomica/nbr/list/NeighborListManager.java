package etomica.nbr.list;

import java.io.Serializable;

import etomica.action.Action;
import etomica.action.AtomAction;
import etomica.action.BoxImposePbc;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IIntegratorNonintervalListener;
import etomica.api.IPotential;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomIteratorTreeBox;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.cell.ApiAACell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.PotentialArray;
import etomica.space.Space;
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
        Action, AgentSource, Serializable {

    /**
     * Configures instance for use by the given PotentialMaster.
     */
    public NeighborListManager(PotentialMasterList potentialMasterList, double range, 
            IBox box, Space space) {
        setUpdateInterval(1);
        this.box = box;
        iieCount = updateInterval;
        iterator = new AtomIteratorTreeBox();
        iterator.setBox(box);
        iterator.setDoAllNodes(true);
        neighborCheck = new NeighborCheck(this);
        setPriority(200);
        pbcEnforcer = new BoxImposePbc(space);
        pbcEnforcer.setBox(box);
        pbcEnforcer.setApplyToMolecules(false);
        potentialMaster = potentialMasterList;
        cellNbrIterator = new ApiAACell(potentialMaster.getSpace().D(), range, box);
        agentManager2Body = new AtomAgentManager(this, box);
        agentManager1Body = new AtomAgentManager(new AtomPotential1ListSource(potentialMasterList), box);
        neighborReset = new NeighborReset(this, agentManager2Body, agentManager1Body);
        boxEvent = new BoxEventNeighborsUpdated(box);
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
        iterator.reset();
        for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            int num2Body = 0;
            int num1Body = 0;
            IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials();
            for (int i=0; i<potentials.length; i++) {
                if (potentials[i].nBody() == 1) {
                    num1Body++;
                }
                else {
                    num2Body++;
                }
            }
            
            ((AtomNeighborLists)agentManager2Body.getAgent(atom)).setCapacity(num2Body);
            ((AtomPotentialList)agentManager1Body.getAgent(atom)).setCapacity(num1Body);
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
        pbcEnforcer.actionPerformed();
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

        iterator.allAtoms(neighborCheck);
        if (neighborCheck.needUpdate) {
            if (Debug.ON && Debug.DEBUG_NOW) {
                System.out.println("Updating neighbors");
            }
            if (neighborCheck.unsafe && !quiet) {
                System.err
                        .println("Atoms exceeded the safe neighbor limit");
            }
            pbcEnforcer.actionPerformed();
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

        iterator.allAtoms(neighborReset);
        
        NeighborCellManager cellManager = potentialMaster.getNbrCellManager(box);
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
                    ((AtomNeighborLists)agentManager2Body.getAgent(atom0)).addUpNbr(atom1,i);
                    ((AtomNeighborLists)agentManager2Body.getAgent(atom1)).addDownNbr(atom0,
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
        return ((AtomNeighborLists)agentManager2Body.getAgent(atom)).getUpList();
    }

    public IAtomSet[] getDownList(IAtom atom) {
        return ((AtomNeighborLists)agentManager2Body.getAgent(atom)).getDownList();
    }

    public AtomPotentialList getPotential1BodyList(IAtom atom) {
        return (AtomPotentialList)agentManager1Body.getAgent(atom);
    }
    
    public void dispose() {
        agentManager1Body.dispose();
        agentManager2Body.dispose();
    }

    private static final long serialVersionUID = 1L;
    private int updateInterval;
    private int iieCount;
    private final AtomIteratorTreeBox iterator;
    private final NeighborCheck neighborCheck;
    private final NeighborReset neighborReset;
    private final ApiAACell cellNbrIterator;
    protected final PotentialMasterList potentialMaster;
    private int priority;
    private BoxImposePbc pbcEnforcer;
    private boolean quiet;
    private final AtomAgentManager agentManager2Body;
    private final AtomAgentManager agentManager1Body;
    protected IBox box;
    private NeighborCriterion[] oldCriteria;
    protected final BoxEventNeighborsUpdated boxEvent;

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

        public NeighborReset(NeighborListManager manager, AtomAgentManager agentManager2Body,
                AtomAgentManager agentManager1Body) {
            neighborListManager = manager;
            this.agentManager2Body = agentManager2Body;
            this.agentManager1Body = agentManager1Body;
            atomSetSinglet = new AtomSetSinglet();
        }
        
        public void actionPerformed(IAtom atom) {
            final NeighborCriterion[] criterion = neighborListManager.getCriterion(atom.getType());
            ((AtomNeighborLists)agentManager2Body.getAgent(atom)).clearNbrs();
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
                ((AtomPotentialList)agentManager1Body.getAgent(atom)).setIsInteracting(criteria[i].accept(atomSetSinglet),i);
            }
        }
        
        private NeighborListManager neighborListManager;
        private AtomAgentManager agentManager2Body;
        private AtomAgentManager agentManager1Body;
        protected final AtomSetSinglet atomSetSinglet;
    }
    
    public Class getAgentClass() {
        return AtomNeighborLists.class;
    }
    
    public Object makeAgent(IAtom atom) {
        AtomNeighborLists lists = new AtomNeighborLists();
        int num2Body = 0;
        IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials();
        for (int i=0; i<potentials.length; i++) {
            if (potentials[i].nBody() > 1) {
                num2Body++;
            }
        }
        lists.setCapacity(num2Body);
        return lists;
    }
    
    public void releaseAgent(Object agent, IAtom atom) {
        ((AtomNeighborLists)agent).clearNbrs();
    }
    
    public static class AtomPotential1ListSource implements AtomAgentManager.AgentSource, java.io.Serializable {
        private static final long serialVersionUID = 2L;
        protected final PotentialMasterList potentialMaster;

        public AtomPotential1ListSource(PotentialMasterList potentialMaster) {
            this.potentialMaster = potentialMaster;
        }
        
        public Class getAgentClass() {
            return AtomPotentialList.class;
        }
        public void releaseAgent(Object obj, IAtom atom) {}
        public Object makeAgent(IAtom atom) {
            AtomPotentialList lists = new AtomPotentialList();
            int num1Body = 0;
            IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials();
            for (int i=0; i<potentials.length; i++) {
                if (potentials[i].nBody() == 1) {
                    num1Body++;
                }
            }
            
            lists.setCapacity(num1Body);
            return lists;
        }
    }
}
