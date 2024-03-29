/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTestCollective;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.nbr.list.NeighborListManager;

public class AtomTestLiquidAtomic implements AtomTestCollective, AtomLeafAgentManager.AgentSource<Boolean> {

    public AtomTestLiquidAtomic(NeighborListManager neighborManager, Box box) {
        leafList = box.getLeafList();
        nbrListManager = neighborManager;
        setMaxNbrsVapor(80);
        agentManager = new AtomLeafAgentManager<Boolean>(this, box);
    }

    public void setMaxNbrsVapor(int newMaxNbrsVapor) {
        maxNbrsVapor = newMaxNbrsVapor;
    }

    public int getMaxNbrsVapor() {
        return maxNbrsVapor;
    }

    public void resetTest() {
        //color all atoms according to their type
        int nLeaf = leafList.size();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtom atom = leafList.get(iLeaf);
            int nbrs = nbrListManager.numAtomNbrsUp[iLeaf] +
                    nbrListManager.numAtomNbrsDn[iLeaf];
            agentManager.setAgent(atom, nbrs > maxNbrsVapor);
        }
    }

    public boolean test(IAtom a) {
        Boolean b = agentManager.getAgent(a);
        return b == null ? false : b;
    }

    public Boolean makeAgent(IAtom a, Box agentBox) {
        return null;
    }

    public void releaseAgent(Boolean agent, IAtom atom, Box agentBox) {
    }

    private final NeighborListManager nbrListManager;
    private final IAtomList leafList;
    protected int maxNbrsVapor;
    protected AtomLeafAgentManager<Boolean> agentManager;
}
