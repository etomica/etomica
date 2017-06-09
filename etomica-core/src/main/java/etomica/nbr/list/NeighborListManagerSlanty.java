/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.PotentialMasterList.NeighborListAgentSource;
import etomica.potential.IPotential;
import etomica.potential.PotentialArray;
import etomica.space.Space;

/**
 * Implements neighbor listing for a slanty box.  Because using slanty cells
 * would be hard, we avoid them entirely.  Neighbors are found by looping over
 * all pairs of atoms.  This is slow, but you're running a solid, right?  So you
 * don't really mind.
 * 
 * @author Andrew Schultz
 */
public class NeighborListManagerSlanty extends NeighborListManager {

    /**
     * Configures instance for use by the given PotentialMaster.
     */
    public NeighborListManagerSlanty(PotentialMasterList potentialMasterList, double range,
                                     Box box, Space space) {
        super(potentialMasterList, range, box, space);
        pair = new AtomPair();
    }

    /**
     * Reassigns all interacting atoms to cells, then loops over all atom
     * pairs, determines for each pair whether a potential applies to it,
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
        
        IAtomList atomList = box.getLeafList();
        for (int iAtom=0; iAtom<atomList.getAtomCount()-1; iAtom++) {
            IAtom atom0 = atomList.getAtom(iAtom);
            pair.atom0 = atom0;
            PotentialArray potentialArray = potentialMaster.getRangedPotentials(atom0.getType());
            IPotential[] potentials = potentialArray.getPotentials();
            NeighborCriterion[] criteria = potentialArray.getCriteria();

            for (int jAtom=iAtom+1; jAtom<atomList.getAtomCount(); jAtom++) {
        
                IAtom atom1 = atomList.getAtom(jAtom);
                pair.atom1 = atom1;
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
        if (agentManager2Body.getAgent(atom) == null) {
            // we're getting called before our own makeAgent (we have no
            // control over order here).  make the agent now and then use it
            // again from makeAgent (ugh).  this depends on AtomAgentManager
            // nulling out agents for removed atoms.
            agentManager2Body.setAgent(atom, makeAgent(atom, box));
        }
        pair.atom0 = atom;
        IAtomList atomList = box.getLeafList();
        PotentialArray potentialArray = potentialMaster.getRangedPotentials(atom.getType());
        IPotential[] potentials = potentialArray.getPotentials();
        NeighborCriterion[] criteria = potentialArray.getCriteria();
        for (int jAtom=0; jAtom<atomList.getAtomCount(); jAtom++) {
            if (jAtom == atom.getLeafIndex()) {
                continue;
            }
            IAtom atom1 = atomList.getAtom(jAtom);
            if (jAtom < atom.getLeafIndex()) {
                pair.atom1 = atom;
                pair.atom0 = atom1;
            }
            else {
                pair.atom0 = atom;
                pair.atom1 = atom1;
            }
            for (int i = 0; i < potentials.length; i++) {
                if (potentials[i].nBody() < 2) {
                    continue;
                }
                if (criteria[i].accept(pair)) {
                    agentManager2Body.getAgent(pair.atom0).addUpNbr(pair.atom1,i);
                    agentManager2Body.getAgent(pair.atom1).addDownNbr(pair.atom0,
                            potentialMaster.getRangedPotentials(atom1.getType()).getPotentialIndex(potentials[i]));
                }
            }
        }
    }
    private static final long serialVersionUID = 1L;
    protected final AtomPair pair;

    /**
     * Constructs instances of NeighborListManagerSlanty on behalf of the
     * PotentialMaster
     */
    public static class NeighborListSlantyAgentSource extends NeighborListAgentSource {
        public NeighborListSlantyAgentSource(double range, Space space) {
            super(range, space);
        }
        
        public NeighborListManager makeAgent(Box box) {
            return new NeighborListManagerSlanty(potentialMaster, range, box, space);
        }
    }
    
}
