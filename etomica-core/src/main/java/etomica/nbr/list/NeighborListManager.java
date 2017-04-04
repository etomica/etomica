/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import java.io.Serializable;

import etomica.action.BoxImposePbc;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IPotential;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomSetSinglet;
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
public class NeighborListManager implements IIntegratorListener, AgentSource<AtomNeighborLists> {

    /**
     * Configures instance for use by the given PotentialMaster.
     */
    public NeighborListManager(PotentialMasterList potentialMasterList, double range, 
            IBox box, ISpace space) {
        setUpdateInterval(1);
        this.box = box;
        iieCount = updateInterval;
        pbcEnforcer = new BoxImposePbc(space);
        pbcEnforcer.setBox(box);
        pbcEnforcer.setApplyToMolecules(false);
        potentialMaster = potentialMasterList;
        cellNbrIterator = new ApiAACell(space.D(), range, box);
        cell1ANbrIterator = new Api1ACell(space.D(), range, potentialMasterList.getCellAgentManager());
        agentManager2Body = new AtomLeafAgentManager<AtomNeighborLists>(this, box,AtomNeighborLists.class);
        AtomPotential1ListSource source1 = new AtomPotential1ListSource(potentialMasterList);
        agentManager1Body = new AtomLeafAgentManager<AtomPotentialList>(source1, box, AtomPotentialList.class);
        source1.setAgentManager(agentManager1Body);
        initialized = false;
        doApplyPBC = true;
        atomSetSinglet = new AtomSetSinglet();
        eventManager = new NeighborListEventManager();
    }

    public void setDoApplyPBC(boolean newDoApplyPBC) {
        doApplyPBC = newDoApplyPBC;
    }

    public boolean getDoApplyPBC() {
        return doApplyPBC;
    }

    public void updateLists() {
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int j=0; j<nLeaf; j++) {
            IAtom atom = leafList.getAtom(j);
            IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials();

            agentManager2Body.getAgent(atom).setCapacity(potentials.length);
            agentManager1Body.getAgent(atom).setCapacity(potentials.length);
        }
    }


    public void integratorInitialized(IIntegratorEvent e) {
        reset();
    }

    public void integratorStepFinished(IIntegratorEvent e) {
        if (--iieCount == 0) {
            updateNbrsIfNeeded();
            iieCount = updateInterval;
        }
    }

    public void integratorStepStarted(IIntegratorEvent e) {}

    /**
     * For each box in the array, applies central image, 
     * resets neighbors of all atoms, and sets up all neighbor
     * lists.
     */
    public void reset() {
        // the NeighborCellManager might not have existed during construction
        // so we couldn't set the lattice.  It better exist now.
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

        boolean needUpdate = false;
        boolean unsafe = false;
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int j=0; j<nLeaf; j++) {
            IAtom atom = leafList.getAtom(j);
            final NeighborCriterion[] criterion = potentialMaster.getRangedPotentials(atom.getType()).getCriteria();
            for (int i = 0; i < criterion.length; i++) {
                if (criterion[i].needUpdate(atom)) {
                    needUpdate = true;
                    if (quiet && (!Debug.ON || !Debug.DEBUG_NOW)) {
                        break;
                    }
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
        
        if (needUpdate) {
            if (Debug.ON && Debug.DEBUG_NOW) {
                System.out.println("Updating neighbors");
            }
            if (unsafe) {
                numUnsafe++;
                if (numUnsafe == 1 || (Long.toString(numUnsafe).matches("10*"))) {
                    System.err.print("Atoms exceeded the safe neighbor limit");
                    if (numUnsafe > 1) {
                        System.err.print(" ("+numUnsafe+" times)");
                    }
                    System.err.println();
                }
            }
            if (doApplyPBC) {
                pbcEnforcer.actionPerformed();
            }
            neighborSetup();
            numUpdates++;
            eventManager.neighborsUpdated();
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

    /**
     * Returns the number of times the neighbor lists have been updated.
     */
    public int getNumUpdates() {
        return numUpdates;
    }

    public NeighborCriterion[] getCriterion(IAtomType atomType) {
        return potentialMaster.getRangedPotentials(atomType).getCriteria();
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
     */
    protected void neighborSetup() {

        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        // reset criteria
        for (int j=0; j<nLeaf; j++) {
            IAtom atom = leafList.getAtom(j);
            final NeighborCriterion[] criterion = getCriterion(atom.getType());
            agentManager2Body.getAgent(atom).clearNbrs();
            for (int i = 0; i < criterion.length; i++) {
                criterion[i].reset(atom);
            }

            PotentialArray potentialArray = potentialMaster.getRangedPotentials(atom.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            NeighborCriterion[] criteria = potentialArray.getCriteria();

            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() != 1) {
                    continue;
                }
                atomSetSinglet.atom = atom;
                agentManager1Body.getAgent(atom).setIsInteracting(criteria[i].accept(atomSetSinglet),i);
            }
        }
        
        NeighborCellManager cellManager = potentialMaster.getNbrCellManager(box);
        cellManager.setDoApplyPBC(!doApplyPBC);
        cellManager.assignCellAll();

        cellNbrIterator.reset();
        //TODO change looping scheme so getPotentials isn't called for every pair
        //consider doing this by introducing ApiNested interface, with hasNextInner and hasNextOuter methods
        for (IAtomList pair = cellNbrIterator.nextPair(); pair != null;
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
                    agentManager2Body.getAgent(atom0).addUpNbr(atom1,i);
                    agentManager2Body.getAgent(atom1).addDownNbr(atom0,
                            potentialMaster.getRangedPotentials(atom1.getType()).getPotentialIndex(potentials[i]));
                }
            }
        }
        initialized = true;
    }

    /**
     * Constructs neighbor lists for the given atom
     */
    public void addAtomNotify(IAtom atom) {
        if (!initialized) {
            // the simulation hasn't started yet.  just wait for neighborSetup
            // to get called.  It can do everything at once and can be sure
            // things are set up properly (like calling setBox on the criteria).
            return;
        }
        
        PotentialArray potentialArray = potentialMaster.getRangedPotentials(atom.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        NeighborCriterion[] criteria = potentialArray.getCriteria();

        if (agentManager1Body.getAgent(atom) == null) {
        	agentManager1Body.setAgent(atom, new AtomPotentialList());
        }
        for (int i = 0; i < potentials.length; i++) {
            if (potentials[i].nBody() != 1) {
                continue;
            }
            atomSetSinglet.atom = atom;
            agentManager1Body.getAgent(atom).setIsInteracting(criteria[i].accept(atomSetSinglet),i);
        }
        
        if (agentManager2Body.getAgent(atom) == null) {
            // we're getting called before our own makeAgent (we have no
            // control over order here).  make the agent now and then use it
            // again from makeAgent (ugh).  this depends on AtomAgentManager
            // nulling out agents for removed atoms.
            agentManager2Body.setAgent(atom, makeAgent(atom, box));
        }
        cell1ANbrIterator.setBox(box);
        cell1ANbrIterator.setTarget(atom);
        cell1ANbrIterator.reset();
        for (IAtomList pair = cell1ANbrIterator.next(); pair != null;
             pair = cell1ANbrIterator.next()) {
            IAtom atom1 = pair.getAtom(1);
            if (atom1 == atom) atom1 = pair.getAtom(0);
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() < 2) {
                    continue;
                }
                if (criteria[i].accept(pair)) {
                    agentManager2Body.getAgent(atom).addUpNbr(atom1,i);
                    agentManager2Body.getAgent(atom1).addDownNbr(atom,
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
    
    public IAtomList[] getUpList(IAtom atom) {
        return agentManager2Body.getAgent(atom).getUpList();
    }

    public IAtomList[] getDownList(IAtom atom) {
        return agentManager2Body.getAgent(atom).getDownList();
    }

    public AtomPotentialList getPotential1BodyList(IAtom atom) {
        return agentManager1Body.getAgent(atom);
    }
    
    public void dispose() {
        agentManager1Body.dispose();
        agentManager2Body.dispose();
    }

    public NeighborListEventManager getEventManager() {
        return eventManager;
    }
    
    private static final long serialVersionUID = 1L;
    private int updateInterval;
    private int iieCount;
    protected final AtomSetSinglet atomSetSinglet;
    private final ApiAACell cellNbrIterator;
    private final Api1ACell cell1ANbrIterator;
    protected final PotentialMasterList potentialMaster;
    private BoxImposePbc pbcEnforcer;
    private boolean quiet;
    protected long numUnsafe;
    protected final AtomLeafAgentManager<AtomNeighborLists> agentManager2Body;
    protected final AtomLeafAgentManager<AtomPotentialList> agentManager1Body;
    private NeighborListEventManager eventManager;
    protected IBox box;
    private NeighborCriterion[] oldCriteria;
    protected boolean initialized;
    protected boolean doApplyPBC;
    protected int numUpdates;

    public AtomNeighborLists makeAgent(IAtom atom, IBox agentBox) {
        if (initialized) {
            AtomNeighborLists oldAgent = agentManager2Body.getAgent(atom);
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
    
    public void releaseAgent(AtomNeighborLists agent, IAtom atom, IBox agentBox) {
        // we need to remove this atom from the neighbor lists of its neighbors.
        AtomNeighborLists nbrLists = agent;
        IAtomList[] upDnLists = nbrLists.getUpList();
        for (int i=0; i<upDnLists.length; i++) {
            int nNbrs = upDnLists[i].getAtomCount();
            for (int j=0; j<nNbrs; j++) {
                IAtom jAtom = upDnLists[i].getAtom(j);
                AtomNeighborLists jNbrLists = agentManager2Body.getAgent(jAtom);
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
                IAtom jAtom = upDnLists[i].getAtom(j);
                AtomNeighborLists jNbrLists = agentManager2Body.getAgent(jAtom);
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

    public static class AtomPotential1ListSource implements AtomLeafAgentManager.AgentSource<AtomPotentialList> {
        protected final PotentialMasterList potentialMaster;
        protected AtomLeafAgentManager<AtomPotentialList> manager;

        public AtomPotential1ListSource(PotentialMasterList potentialMaster) {
            this.potentialMaster = potentialMaster;
        }
        
        public void setAgentManager(AtomLeafAgentManager<AtomPotentialList> manager) {
        	this.manager = manager;
        }

        public void releaseAgent(AtomPotentialList obj, IAtom atom, IBox agentBox) {}
        public AtomPotentialList makeAgent(IAtom atom, IBox agentBox) {
        	if (manager != null) {
	            AtomPotentialList oldAgent = manager.getAgent(atom);
	            if (oldAgent != null) {
	                // NeighborCellManager got notified first and we already made the
	                // agent (and found the neighbors!).  Return that now.
	                return oldAgent;
	            }
        	}
            AtomPotentialList lists = new AtomPotentialList();
            IPotential[] potentials = potentialMaster.getRangedPotentials(atom.getType()).getPotentials();
            
            lists.setCapacity(potentials.length);
            return lists;
        }
    }
}
