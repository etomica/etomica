/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.AtomFilterCollective;
import etomica.atom.AtomLeafAgentManager;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;

public class AtomFilterLiquidAtomic implements AtomFilterCollective, AtomLeafAgentManager.AgentSource<Boolean> {
    
    public AtomFilterLiquidAtomic(PotentialMasterList potentialMaster, Box box) {
        leafList = box.getLeafList();
        nbrListManager = potentialMaster.getNeighborManager(box);
        setMaxNbrsVapor(80);
        agentManager = new AtomLeafAgentManager<Boolean>(this, box, Boolean.class);
    }

    public void setMaxNbrsVapor(int newMaxNbrsVapor) {
        maxNbrsVapor = newMaxNbrsVapor;
    }
    
    public int getMaxNbrsVapor() {
        return maxNbrsVapor;
    }
    
    public void resetFilter() {
		//color all atoms according to their type
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom atom = leafList.getAtom(iLeaf);
            int nbrs = nbrListManager.getUpList(atom)[0].getAtomCount() +
                       nbrListManager.getDownList(atom)[0].getAtomCount();
            agentManager.setAgent(atom, nbrs > maxNbrsVapor);
        }
    }
    
    public boolean accept(IAtom a) {
        Boolean b = agentManager.getAgent(a);
        return b == null ? false : b;
    }

    public boolean accept(IMolecule mole) {
        return false;
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
