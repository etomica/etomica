/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.action.BoxImposePbc;
import etomica.atom.*;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.ApiAACell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.IPotentialAtomic;
import etomica.space.Space;
import etomica.util.Debug;

import java.util.List;

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
public class NeighborListManager implements IntegratorListener, AgentSource<AtomNeighborLists> {

    private static final long serialVersionUID = 1L;
    protected final AtomSetSinglet atomSetSinglet;
    protected final PotentialMasterList potentialMaster;
    protected final AtomLeafAgentManager<AtomNeighborLists> agentManager2Body;
    protected final AtomLeafAgentManager<AtomPotentialList> agentManager1Body;
    private final ApiAACell cellNbrIterator;
    private final Api1ACell cell1ANbrIterator;
    protected long numUnsafe;
    protected Box box;
    protected boolean initialized;
    protected boolean doApplyPBC;
    protected int numUpdates;
    private int updateInterval;
    private int iieCount;
    private BoxImposePbc pbcEnforcer;
    private boolean quiet;
    private NeighborListEventManager eventManager;
    private NeighborCriterion[] oldCriteria;
    
    /**
     * Configures instance for use by the given PotentialMaster.
     */
    public NeighborListManager(PotentialMasterList potentialMasterList, double range, Box box) {
        setUpdateInterval(1);
        Space space = box.getSpace();
        this.box = box;
        iieCount = updateInterval;
        pbcEnforcer = new BoxImposePbc(space);
        pbcEnforcer.setBox(box);
        pbcEnforcer.setApplyToMolecules(false);
        potentialMaster = potentialMasterList;
        cellNbrIterator = new ApiAACell(space.D(), range, box);
        cell1ANbrIterator = new Api1ACell(space.D(), range, potentialMasterList.getCellAgentManager());
        agentManager2Body = new AtomLeafAgentManager<AtomNeighborLists>(this, box);
        AtomPotential1ListSource source1 = new AtomPotential1ListSource(potentialMasterList);
        agentManager1Body = new AtomLeafAgentManager<AtomPotentialList>(source1, box);
        source1.setAgentManager(agentManager1Body);
        initialized = false;
        doApplyPBC = true;
        atomSetSinglet = new AtomSetSinglet();
        eventManager = new NeighborListEventManager();
    }

    public boolean getDoApplyPBC() {
        return doApplyPBC;
    }

    public void setDoApplyPBC(boolean newDoApplyPBC) {
        doApplyPBC = newDoApplyPBC;
    }
    
    public void integratorInitialized(IntegratorEvent e) {
        reset();
    }

    public void integratorStepFinished(IntegratorEvent e) {
        if (--iieCount == 0) {
            updateNbrsIfNeeded();
            iieCount = updateInterval;
        }
    }

    public void integratorStepStarted(IntegratorEvent e) {}

    /**
     * For each box in the array, applies central image,
     * resets neighbors of all atoms, and sets up all neighbor
     * lists.
     */
    public void reset() {
        // the NeighborCellManager might not have existed during construction
        // so we couldn't set the lattice.  It better exist now.
        cellNbrIterator.setLattice(potentialMaster.getNbrCellManager(box).getLattice());

        potentialMaster.setBoxForCriteria(box);
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

        potentialMaster.setBoxForCriteria(box);

        boolean needUpdate = false;
        boolean unsafe = false;
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.size();
        for (int j=0; j<nLeaf; j++) {
            IAtom atom = leafList.get(j);

            final NeighborCriterion[] criterion = potentialMaster.getCriteria(atom.getType());
            for (int i = 0; i < criterion.length; i++) {
                if (criterion[i] == null || !criterion[i].needUpdate(atom)) continue;
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

    public NeighborCriterion[] getCriteria(AtomType atomType) {
        return potentialMaster.getCriteria(atomType);
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
        // reset criteria
        for (IAtom atom : leafList) {
            final NeighborCriterion[] criterion = potentialMaster.getCriteria(atom.getType());
            agentManager2Body.getAgent(atom).clearNbrs();
            for (NeighborCriterion aCriterion : criterion) {
                if (aCriterion != null) aCriterion.reset(atom);
            }

            List<NeighborCriterion> criteria = potentialMaster.getCriteria1Body(atom.getType());

            for (int i = 0; i < criteria.size(); i++) {
                atomSetSinglet.atom = atom;
                agentManager1Body.getAgent(atom).setIsInteracting(criteria.get(i).accept(atomSetSinglet), i);
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
            IAtom atom0 = pair.get(0);
            IAtom atom1 = pair.get(1);
            int type1 = atom1.getType().getIndex();
            NeighborCriterion criteria = potentialMaster.getCriteria(atom0.getType())[type1];
            if (criteria != null && criteria.accept(pair)) {
                agentManager2Body.getAgent(atom0).addUpNbr(atom1, type1);
                agentManager2Body.getAgent(atom1).addDownNbr(atom0, atom0.getType().getIndex());
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

        List<NeighborCriterion> criteria1 = potentialMaster.getCriteria1Body(atom.getType());

        if (agentManager1Body.getAgent(atom) == null) {
        	agentManager1Body.setAgent(atom, new AtomPotentialList());
        }
        for (int i = 0; i < criteria1.size(); i++) {
            NeighborCriterion c = criteria1.get(i);
            atomSetSinglet.atom = atom;
            agentManager1Body.getAgent(atom).setIsInteracting(c.accept(atomSetSinglet), i);
        }

        NeighborCriterion[] criteria = potentialMaster.getCriteria(atom.getType());

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
            IAtom atom1 = pair.get(1);
            if (atom1 == atom) atom1 = pair.get(0);
            NeighborCriterion criterion = criteria[atom1.getType().getIndex()];
            if (criterion != null && criterion.accept(pair)) {
                agentManager2Body.getAgent(atom).addUpNbr(atom1, atom1.getType().getIndex());
                agentManager2Body.getAgent(atom1).addDownNbr(atom, atom.getType().getIndex());
            }
        }
    }

    public double getRange() {
        return cellNbrIterator.getNbrCellIterator().getNeighborDistance();
    }

    /**
     * Sets the interaction range, which affects the cell-list neighbor iteration
     * used to generate candidate neighbors for neighbor listing.
     */
    public void setRange(double d) {
        cell1ANbrIterator.getNbrCellIterator().setNeighborDistance(d);
        cellNbrIterator.getNbrCellIterator().setNeighborDistance(d);
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

    public AtomNeighborLists makeAgent(IAtom atom, Box agentBox) {
        if (initialized) {
            AtomNeighborLists oldAgent = agentManager2Body.getAgent(atom);
            if (oldAgent != null) {
                // NeighborCellManager got notified first and we already made the
                // agent (and found the neighbors!).  Return that now.
                return oldAgent;
            }
        }
        AtomNeighborLists lists = new AtomNeighborLists();
        IPotentialAtomic[] potentials = potentialMaster.getRangedPotentials(atom.getType());
        lists.setCapacity(potentials.length);
        return lists;
    }
    
    public void releaseAgent(AtomNeighborLists agent, IAtom atom, Box agentBox) {
        // we need to remove this atom from the neighbor lists of its neighbors.
        AtomNeighborLists nbrLists = agent;
        IAtomList[] upDnLists = nbrLists.getUpList();
        for (int i=0; i<upDnLists.length; i++) {
            int nNbrs = upDnLists[i].size();
            for (int j=0; j<nNbrs; j++) {
                IAtom jAtom = upDnLists[i].get(j);
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
            int nNbrs = upDnLists[i].size();
            for (int j=0; j<nNbrs; j++) {
                IAtom jAtom = upDnLists[i].get(j);
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

        public void releaseAgent(AtomPotentialList obj, IAtom atom, Box agentBox) {}
        public AtomPotentialList makeAgent(IAtom atom, Box agentBox) {
        	if (manager != null) {
	            AtomPotentialList oldAgent = manager.getAgent(atom);
	            if (oldAgent != null) {
	                // NeighborCellManager got notified first and we already made the
	                // agent (and found the neighbors!).  Return that now.
	                return oldAgent;
	            }
        	}
            AtomPotentialList lists = new AtomPotentialList();
            List<IPotentialAtomic> potentials = potentialMaster.getRangedPotentials1Body(atom.getType());
            lists.setCapacity(potentials.size());
            return lists;
        }
    }
}
