/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.list;

import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.nbr.NeighborCriterion;

import java.util.List;

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
                                     Box box) {
        super(potentialMasterList, range, box);
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
        int nLeaf = leafList.size();
        // reset criteria
        for (int j=0; j<nLeaf; j++) {
            IAtom atom = leafList.get(j);
            final NeighborCriterion[] criterion = potentialMaster.getCriteria(atom.getType());
            agentManager2Body.getAgent(atom).clearNbrs();
            for (int i = 0; i < criterion.length; i++) {
                criterion[i].reset(atom);
            }

            List<NeighborCriterion> criteria = potentialMaster.getCriteria1Body(atom.getType());

            for (int i = 0; i < criteria.size(); i++) {
                agentManager1Body.getAgent(atom).setIsInteracting(criteria.get(i).accept(atom, null), i);
            }
        }
        
        IAtomList atomList = box.getLeafList();
        for (int iAtom = 0; iAtom<atomList.size()-1; iAtom++) {
            IAtom atom0 = atomList.get(iAtom);
            pair.atom0 = atom0;
            NeighborCriterion[] criteria = potentialMaster.getCriteria(atom0.getType());

            for (int jAtom = iAtom+1; jAtom<atomList.size(); jAtom++) {
        
                IAtom atom1 = atomList.get(jAtom);
                NeighborCriterion c = criteria[atom1.getType().getIndex()];
                if (c == null) continue;
                pair.atom1 = atom1;
                if (c.accept(atom0, atom1)) {
                    agentManager2Body.getAgent(atom0).addUpNbr(atom1, atom0.getType().getIndex());
                    agentManager2Body.getAgent(atom1).addDownNbr(atom0, atom1.getType().getIndex());
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
        NeighborCriterion[] criteria = potentialMaster.getCriteria(atom.getType());
        IAtom firstAtom, secondAtom;
        for (int jAtom = 0; jAtom<atomList.size(); jAtom++) {
            if (jAtom == atom.getLeafIndex()) {
                continue;
            }
            IAtom atom1 = atomList.get(jAtom);
            NeighborCriterion c = criteria[atom1.getType().getIndex()];
            if (c == null) continue;
            if (jAtom < atom.getLeafIndex()) {
                firstAtom = atom;
                secondAtom = atom1;
            }
            else {
                firstAtom = atom;
                secondAtom = atom1;
            }

            if (c.accept(firstAtom, secondAtom)) {
                agentManager2Body.getAgent(firstAtom).addUpNbr(secondAtom, pair.atom0.getType().getIndex());
                agentManager2Body.getAgent(secondAtom).addDownNbr(firstAtom, pair.atom1.getType().getIndex());
            }
        }
    }
    protected final AtomPair pair;

    /**
     * Constructs instances of NeighborListManagerSlanty on behalf of the
     * PotentialMaster
     */
    public static class NeighborListSlantyAgentSource extends NeighborListAgentSource {
        public NeighborListSlantyAgentSource(double range) {
            super(range);
        }
        
        public NeighborListManager makeAgent(Box box) {
            return new NeighborListManagerSlanty(potentialMaster, range, box);
        }
    }
    
}
