package etomica.modules.droplet;

import java.util.Random;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.atom.AtomFilterCollective;
import etomica.atom.AtomLeafAgentManager;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;

public class AtomFilterLiquidAtomic implements AtomFilterCollective, AtomLeafAgentManager.AgentSource {
    
    public AtomFilterLiquidAtomic(PotentialMasterList potentialMaster, IBox box) {
        leafList = box.getLeafList();
        nbrListManager = potentialMaster.getNeighborManager(box);
        setMaxNbrsVapor(80);
        agentManager = new AtomLeafAgentManager(this, box);
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
        int nLiquid = 0;
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom atom = leafList.getAtom(iLeaf);
            int nbrs = nbrListManager.getUpList(atom)[0].getAtomCount() +
                       nbrListManager.getDownList(atom)[0].getAtomCount();
            nLiquid += nbrs > maxNbrsVapor ? 1 : 0;
            agentManager.setAgent(atom, nbrs > maxNbrsVapor);
        }
        if (new Random().nextInt(100) == 5) {
            System.out.println("numLiquid "+nLiquid);
        }
    }
    
    public boolean accept(IAtom a) {
        Boolean b = (Boolean)agentManager.getAgent(a);
        return b == null ? false : b;
    }

    public boolean accept(IMolecule mole) {
        return false;
    }

    public Class getAgentClass() {
        return Boolean.class;
    }

    public Object makeAgent(IAtom a) {
        return null;
    }

    public void releaseAgent(Object agent, IAtom atom) {
    }

    private static final long serialVersionUID = 1L;
    private final NeighborListManager nbrListManager;
    private final IAtomList leafList;
    protected int maxNbrsVapor;
    protected AtomLeafAgentManager agentManager;
}