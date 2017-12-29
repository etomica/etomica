/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list.molecule;

import etomica.action.BoxImposePbc;
import etomica.box.Box;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.*;
import etomica.nbr.cell.molecule.Mpi1ACell;
import etomica.nbr.cell.molecule.MpiAACell;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.nbr.molecule.NeighborCriterionMolecular;
import etomica.potential.IPotential;
import etomica.potential.IPotentialMolecular;
import etomica.potential.PotentialArrayMolecular;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.Debug;

import java.io.Serializable;

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
 *
 * @author taitan
 *
 */
public class NeighborListManagerMolecular implements IntegratorListener, MoleculeAgentManager.MoleculeAgentSource, Serializable {

    /**
     * Configures instance for use by the given PotentialMaster.
     */
    public NeighborListManagerMolecular(PotentialMasterListMolecular potentialMasterList, double range,
                                        Box box, Space space) {
        setUpdateInterval(1);
        this.box = box;
        iieCount = updateInterval;
        pbcEnforcer = new BoxImposePbc(space);
        pbcEnforcer.setBox(box);
        pbcEnforcer.setApplyToMolecules(true);
        potentialMaster = potentialMasterList;
        cellNbrIterator = new MpiAACell(space.D(), range, box);
        cell1ANbrIterator = new Mpi1ACell(space.D(), range, potentialMasterList.getCellAgentManager());
        
        agentManager2Body = new MoleculeAgentManager(potentialMasterList.getSimulation(), box, this);
        agentManager1Body = new MoleculeAgentManager(potentialMasterList.getSimulation(), box, new MoleculePotential1ListSource(potentialMasterList));
        initialized = false;
        doApplyPBC = true;
        moleculeSetSinglet = new MoleculeSetSinglet();
        eventManager = new NeighborListEventManagerMolecular();
    }

    public void setDoApplyPBC(boolean newDoApplyPBC) {
        doApplyPBC = newDoApplyPBC;
    }

    public boolean getDoApplyPBC() {
        return doApplyPBC;
    }

    public void updateLists() {
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecule = moleculeList.getMoleculeCount();
        for (int j=0; j<nMolecule; j++) {
            IMolecule molecule = moleculeList.getMolecule(j);
            IPotential[] potentials = potentialMaster.getRangedPotentials(molecule.getType()).getPotentials();

            ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule)).setCapacity(potentials.length);
            ((MoleculePotentialList)agentManager1Body.getAgent(molecule)).setCapacity(potentials.length);
        }
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
     
        NeighborCriterionMolecular[] criteriaArray = potentialMaster.getNeighborCriteria();
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

        NeighborCriterionMolecular[] criteriaArray = potentialMaster.getNeighborCriteria();
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
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecule = moleculeList.getMoleculeCount();
        for (int j=0; j<nMolecule; j++) {
            IMolecule molecule = moleculeList.getMolecule(j);
            final NeighborCriterionMolecular[] criterion = potentialMaster.getRangedPotentials(molecule.getType()).getCriteria();
            for (int i = 0; i < criterion.length; i++) {
                if (criterion[i].needUpdate(molecule)) {
                    needUpdate = true;
                    if (quiet && (!Debug.ON || !Debug.DEBUG_NOW)) {
                        break;
                    }
                    if (criterion[i].unsafe()) {
                        if (Debug.ON && Debug.DEBUG_NOW) {
                            System.out.println("molecule " + molecule
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
                System.err.println("Molecules exceeded the safe neighbor limit");
            }
            if (doApplyPBC) {
                pbcEnforcer.actionPerformed();
            }
            neighborSetup();
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

    public NeighborCriterionMolecular[] getCriterion(ISpecies species) {
        return potentialMaster.getRangedPotentials(species).getCriteria();
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

        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecule = moleculeList.getMoleculeCount();
        // reset criteria
        for (int j=0; j<nMolecule; j++) {
            IMolecule molecule = moleculeList.getMolecule(j);
            final NeighborCriterionMolecular[] criterion = getCriterion(molecule.getType());
            ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule)).clearNbrs();
            for (int i = 0; i < criterion.length; i++) {
                criterion[i].reset(molecule);
            }

            PotentialArrayMolecular potentialArray = potentialMaster.getRangedPotentials(molecule.getType());
            IPotentialMolecular[] potentials = potentialArray.getPotentials();
            NeighborCriterionMolecular[] criteria = potentialArray.getCriteria();
            
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() != 1) {
                    continue;
                }
                moleculeSetSinglet.atom = molecule;
                ((MoleculePotentialList)agentManager1Body.getAgent(molecule)).setIsInteracting(criteria[i].accept(moleculeSetSinglet),i);
            }
        }
        
        NeighborCellManagerMolecular cellManager = potentialMaster.getNbrCellManager(box);
        cellManager.setDoApplyPBC(!doApplyPBC);
        cellManager.assignCellAll();
      
        cellNbrIterator.reset();
        //TODO change looping scheme so getPotentials isn't called for every pair
        //consider doing this by introducing ApiNested interface, with hasNextInner and hasNextOuter methods
        for (IMoleculeList pair = cellNbrIterator.nextPair(); pair != null;
             pair = cellNbrIterator.nextPair()) {
            IMolecule molecule0 = pair.getMolecule(0);
            IMolecule molecule1 = pair.getMolecule(1);
            PotentialArrayMolecular potentialArray = potentialMaster.getRangedPotentials(molecule0.getType());
            IPotentialMolecular[] potentials = potentialArray.getPotentials();
            NeighborCriterionMolecular[] criteria = potentialArray.getCriteria();
            
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() < 2) {
                    continue;
                }
                if (criteria[i].accept(pair)) {
                    ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule0)).addUpNbr(molecule1,i);
                    ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule1)).addDownNbr(molecule0,
                            potentialMaster.getRangedPotentials(molecule1.getType()).getPotentialIndex(potentials[i]));
                }
            }
        }
        initialized = true;
    }

    /**
     * Constructs neighbor lists for the given atom
     */
    public void addMoleculeNotify(IMolecule molecule) {
        if (!initialized) {
            // the simulation hasn't started yet.  just wait for neighborSetup
            // to get called.  It can do everything at once and can be sure
            // things are set up properly (like calling setBox on the criteria).
            return;
        }
        if (agentManager2Body.getAgent(molecule) == null) {
            // we're getting called before our own makeAgent (we have no
            // control over order here).  make the agent now and then use it
            // again from makeAgent (ugh).  this depends on AtomAgentManager
            // nulling out agents for removed atoms.
            agentManager2Body.setAgent(molecule, makeAgent(molecule));
        }
        cell1ANbrIterator.setBox(box);
        cell1ANbrIterator.setTarget(molecule);
        cell1ANbrIterator.reset();
        for (IMoleculeList pair = cell1ANbrIterator.next(); pair != null;
             pair = cell1ANbrIterator.next()) {
            IMolecule molecule0 = pair.getMolecule(0);
            IMolecule molecule1 = pair.getMolecule(1);
            PotentialArrayMolecular potentialArray = potentialMaster.getRangedPotentials(molecule0.getType());
            IPotentialMolecular[] potentials = potentialArray.getPotentials();
            NeighborCriterionMolecular[] criteria = potentialArray.getCriteria();
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() < 2) {
                    continue;
                }
                if (criteria[i].accept(pair)) {
                    ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule0)).addUpNbr(molecule1,i);
                    ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule1)).addDownNbr(molecule0,
                            potentialMaster.getRangedPotentials(molecule1.getType()).getPotentialIndex(potentials[i]));
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
    
    public IMoleculeList[] getUpList(IMolecule molecule) {
        return ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule)).getUpList();
    }

    public IMoleculeList[] getDownList(IMolecule molecule) {
        return ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule)).getDownList();
    }

    public MoleculePotentialList getPotential1BodyList(IMolecule molecule) {
        return (MoleculePotentialList)agentManager1Body.getAgent(molecule);
    }
    
    public void dispose() {
        agentManager1Body.dispose();
        agentManager2Body.dispose();
    }

    public NeighborListEventManagerMolecular getEventManager() {
        return eventManager;
    }
    
    private static final long serialVersionUID = 1L;
    private int updateInterval;
    private int iieCount;
    protected final MoleculeSetSinglet moleculeSetSinglet;
    private final MpiAACell cellNbrIterator;
    private final Mpi1ACell cell1ANbrIterator;
    protected final PotentialMasterListMolecular potentialMaster;
    private BoxImposePbc pbcEnforcer;
    private boolean quiet;
    protected final MoleculeAgentManager agentManager2Body;
    protected final MoleculeAgentManager agentManager1Body;
    private NeighborListEventManagerMolecular eventManager;
    protected Box box;
    private NeighborCriterionMolecular[] oldCriteria;
    protected boolean initialized;
    protected boolean doApplyPBC;

    public Object makeAgent(IMolecule molecule) {
        if (initialized) {
            Object oldAgent = agentManager2Body.getAgent(molecule);
            if (oldAgent != null) {
                // NeighborCellManager got notified first and we already made the
                // agent (and found the neighbors!).  Return that now.
                return oldAgent;
            }
        }
        MoleculeNeighborLists lists = new MoleculeNeighborLists();
        IPotential[] potentials = potentialMaster.getRangedPotentials(molecule.getType()).getPotentials();
        lists.setCapacity(potentials.length);
        return lists;
    }
    
    public void releaseAgent(Object agent, IMolecule molecule) {
        // we need to remove this atom from the neighbor lists of its neighbors.
        MoleculeNeighborLists nbrLists = (MoleculeNeighborLists)agent;
        IMoleculeList[] upDnLists = nbrLists.getUpList();
        for (int i=0; i<upDnLists.length; i++) {
            int nNbrs = upDnLists[i].getMoleculeCount();
            for (int j=0; j<nNbrs; j++) {
                IMolecule jMolecule = upDnLists[i].getMolecule(j);
                MoleculeNeighborLists jNbrLists = (MoleculeNeighborLists)agentManager2Body.getAgent(jMolecule);
                MoleculeArrayList[] jDnLists = jNbrLists.downList;
                for (int k=0; k<jDnLists.length; k++) {
                    int idx = jDnLists[k].indexOf(molecule);
                    if (idx > -1) {
                        jDnLists[k].removeAndReplace(idx);
                    }
                }
            }
        }
        upDnLists = nbrLists.getDownList();
        for (int i=0; i<upDnLists.length; i++) {
            int nNbrs = upDnLists[i].getMoleculeCount();
            for (int j=0; j<nNbrs; j++) {
                IMolecule jMolecule = upDnLists[i].getMolecule(j);
                MoleculeNeighborLists jNbrLists = (MoleculeNeighborLists)agentManager2Body.getAgent(jMolecule);
                MoleculeArrayList[] jUpLists = jNbrLists.upList;
                for (int k=0; k<jUpLists.length; k++) {
                    int idx = jUpLists[k].indexOf(molecule);
                    if (idx > -1) {
                        jUpLists[k].removeAndReplace(idx);
                    }
                }
            }
        }
        nbrLists.clearNbrs();
    }

    public static class MoleculePotential1ListSource implements MoleculeAgentManager.MoleculeAgentSource, java.io.Serializable {
        private static final long serialVersionUID = 2L;
        protected final PotentialMasterListMolecular potentialMaster;

        public MoleculePotential1ListSource(PotentialMasterListMolecular potentialMaster) {
            this.potentialMaster = potentialMaster;
        }
        
        public void releaseAgent(Object obj, IMolecule molecule) {}
        public Object makeAgent(IMolecule molecule) {
            MoleculePotentialList lists = new MoleculePotentialList();
            IPotential[] potentials = potentialMaster.getRangedPotentials(molecule.getType()).getPotentials();
            
            lists.setCapacity(potentials.length);
            return lists;
        }

    }


}
