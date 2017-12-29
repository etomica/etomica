/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.AtomFilterCollective;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;

/**
 */
public class AtomFilterChainLength implements AtomFilterCollective, AtomLeafAgentManager.AgentSource<AtomFilterChainLength.LengthAgent> {

    public AtomFilterChainLength(AtomLeafAgentManager<IAtom[]> aam) {
        agentManager = aam;
    }

    public LengthAgent makeAgent(IAtom a, Box agentBox) {
        return new LengthAgent();
    }

    public void releaseAgent(LengthAgent agent, IAtom atom, Box agentBox) {}

    public void resetFilter() {
        maxChainLength = 0;
        // untag all the Atoms
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            chainLengthManager.getAgent(leafList.getAtom(i)).chainLength = 0;
        }

        for (int i=0; i<nLeaf; i++) {
            IAtom a = leafList.getAtom(i);
            // if an Atom has a chain length, it was already counted as part of 
            // another chain
            if (chainLengthManager.getAgent(a).chainLength > 0) continue;

            int chainLength = recursiveTag(a, -1);
            // re-tag with the actual chain length
            recursiveTag(a, chainLength);
            if (chainLength > maxChainLength) {
                maxChainLength = chainLength;
            }
        }
    }

    protected int recursiveTag(IAtom a, int chainLength) {
        chainLengthManager.getAgent(a).chainLength = chainLength;

        IAtom[] nbrs = agentManager.getAgent(a);

        int ctr = 1;
        
        // count all the bonded partners
        for(int i=0; i<nbrs.length; i++) {
            if(nbrs[i] == null) continue;
            if(chainLengthManager.getAgent(nbrs[i]).chainLength == chainLength) {
                // this Atom was already counted as being within this chain
                // so skip it
                continue;
            }
            // count this Atom and all of its bonded partners
            ctr += recursiveTag(nbrs[i], chainLength);
        }
        return ctr;
    }

    public Box getBox() {
        return box;
    }

    public void setBox(Box box) {
        this.box = box;
        if (chainLengthManager != null) {
            // allow old agentManager to de-register itself as a BoxListener
            chainLengthManager.dispose();
        }
        chainLengthManager = new AtomLeafAgentManager<LengthAgent>(this,box);
    }

    public boolean test(IAtom a) {
        return chainLengthManager.getAgent(a).chainLength == maxChainLength;
    }

    protected Box box;
    protected AtomLeafAgentManager<LengthAgent> chainLengthManager;
    protected AtomLeafAgentManager<IAtom[]> agentManager;
    protected int maxChainLength = 0;

    public static class LengthAgent {
        public int chainLength;
    }
}
