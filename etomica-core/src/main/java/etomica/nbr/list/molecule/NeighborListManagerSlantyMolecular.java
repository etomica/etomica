/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list.molecule;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePair;
import etomica.nbr.list.molecule.PotentialMasterListMolecular.NeighborListAgentSourceMolecular;
import etomica.nbr.molecule.NeighborCriterionMolecular;
import etomica.potential.IPotential;
import etomica.potential.IPotentialMolecular;
import etomica.potential.PotentialArrayMolecular;
import etomica.space.Space;

/**
 * Implements neighbor listing for a slanty box.  Because using slanty cells
 * would be hard, we avoid them entirely.  Neighbors are found by looping over
 * all pairs of atoms.  This is slow, but you're running a solid, right?  So you
 * don't really mind.
 * 
 * @author 	Tai Boon Tan
 */
public class NeighborListManagerSlantyMolecular extends NeighborListManagerMolecular {

    /**
     * Configures instance for use by the given PotentialMaster.
     */
    public NeighborListManagerSlantyMolecular(PotentialMasterListMolecular potentialMasterList, double range,
                                              Box box, Space space) {
        super(potentialMasterList, range, box, space);
        pair = new MoleculePair();
    }

    /**
     * Reassigns all interacting atoms to cells, then loops over all atom
     * pairs, determines for each pair whether a potential applies to it,
     * and if so, puts each in the other's neighbor list.
     * Called by updateNbrsIfNeeded, and by reset.
     */
    protected void neighborSetup() {

        IMoleculeList moleculeList = box.getMoleculeList();
        int nLeaf = moleculeList.size();
        // reset criteria
        for (int j=0; j<nLeaf; j++) {
            IMolecule molecule = moleculeList.get(j);
            final NeighborCriterionMolecular[] criterion = getCriterion(molecule.getType());
            ((MoleculeNeighborLists)agentManager2Body.getAgent(molecule)).clearNbrs();
            for (int i = 0; i < criterion.length; i++) {
                criterion[i].reset(molecule);
            }

            PotentialArrayMolecular potentialArray = potentialMaster.getRangedPotentials(molecule.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            NeighborCriterionMolecular[] criteria = potentialArray.getCriteria();

            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() != 1) {
                    continue;
                }
                moleculeSetSinglet.atom = molecule;
                ((MoleculePotentialList)agentManager1Body.getAgent(molecule)).setIsInteracting(criteria[i].accept(moleculeSetSinglet),i);
            }
        }
        
        moleculeList = box.getMoleculeList();
        for (int iMolecule = 0; iMolecule<moleculeList.size()-1; iMolecule++) {
            IMolecule molecule0 = moleculeList.get(iMolecule);
            pair.atom0 = molecule0;
            PotentialArrayMolecular potentialArray = potentialMaster.getRangedPotentials(molecule0.getType());
            IPotentialMolecular[] potentials = potentialArray.getPotentials();
            NeighborCriterionMolecular[] criteria = potentialArray.getCriteria();

            for (int jMolecule = iMolecule+1; jMolecule<moleculeList.size(); jMolecule++) {
        
                IMolecule molecule1 = moleculeList.get(jMolecule);
                pair.atom1 = molecule1;
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
        initialized = true;
    }

    /**
     * Constructs neighbor lists for the given molecule
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
        pair.atom0 = molecule;
        IMoleculeList moleculeList = box.getMoleculeList();
        PotentialArrayMolecular potentialArray = potentialMaster.getRangedPotentials(molecule.getType());
        IPotentialMolecular[] potentials = potentialArray.getPotentials();
        NeighborCriterionMolecular[] criteria = potentialArray.getCriteria();
        for (int jMolecule = 0; jMolecule<moleculeList.size(); jMolecule++) {
            if (jMolecule == molecule.getIndex()) {
                continue;
            }
            IMolecule molecule1 = moleculeList.get(jMolecule);
            if (jMolecule < molecule.getIndex()) {
                pair.atom1 = molecule;
                pair.atom0 = molecule1;
            }
            else {
                pair.atom0 = molecule;
                pair.atom1 = molecule1;
            }
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() < 2) {
                    continue;
                }
                if (criteria[i].accept(pair)) {
                    ((MoleculeNeighborLists)agentManager2Body.getAgent(pair.atom0)).addUpNbr(pair.atom1,i);
                    ((MoleculeNeighborLists)agentManager2Body.getAgent(pair.atom1)).addDownNbr(pair.atom0,
                            potentialMaster.getRangedPotentials(molecule1.getType()).getPotentialIndex(potentials[i]));
                }
            }
        }
    }
    private static final long serialVersionUID = 1L;
    protected final MoleculePair pair;

    /**
     * Constructs instances of NeighborListManagerSlanty on behalf of the
     * PotentialMaster
     */
    public static class NeighborListSlantyAgentSourceMolecular extends NeighborListAgentSourceMolecular {
        public NeighborListSlantyAgentSourceMolecular(double range, Space space) {
            super(range, space);
        }

        public NeighborListManagerMolecular makeAgent(Box box) {
            return new NeighborListManagerSlantyMolecular(potentialMaster, range, box, space);
        }
    }
    
}
