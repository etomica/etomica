package etomica.modules.chainequilibrium;

import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.atom.AtomFilterCollective;
import etomica.atom.AtomLeafAgentManager;

/**
 */
public class AtomFilterChainLength implements AtomFilterCollective, AtomLeafAgentManager.AgentSource {

    public AtomFilterChainLength(AtomLeafAgentManager aam) {
        agentManager = aam;
    }

    public Class getAgentClass() {
        return LengthAgent.class;
    }

    public Object makeAgent(IAtomLeaf a) {
        return new LengthAgent();
    }

    public void releaseAgent(Object agent, IAtomLeaf atom) {}

    public void resetFilter() {
        maxChainLength = 0;
        // untag all the Atoms
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            ((LengthAgent)chainLengthManager.getAgent(leafList.getAtom(i))).chainLength = 0;
        }

        for (int i=0; i<nLeaf; i++) {
            IAtomLeaf a = leafList.getAtom(i);
            // if an Atom has a chain length, it was already counted as part of 
            // another chain
            if (((LengthAgent)chainLengthManager.getAgent(a)).chainLength > 0) continue;

            int chainLength = recursiveTag(a, -1);
            // re-tag with the actual chain length
            recursiveTag(a, chainLength);
            if (chainLength > maxChainLength) {
                maxChainLength = chainLength;
            }
        }
    }

    protected int recursiveTag(IAtomLeaf a, int chainLength) {
        ((LengthAgent)chainLengthManager.getAgent(a)).chainLength = chainLength;

        IAtomLeaf[] nbrs = (IAtomLeaf[])agentManager.getAgent(a);

        int ctr = 1;
        
        // count all the bonded partners
        for(int i=0; i<nbrs.length; i++) {
            if(nbrs[i] == null) continue;
            if(((LengthAgent)chainLengthManager.getAgent(nbrs[i])).chainLength == chainLength) {
                // this Atom was already counted as being within this chain
                // so skip it
                continue;
            }
            // count this Atom and all of its bonded partners
            ctr += recursiveTag(nbrs[i], chainLength);
        }
        return ctr;
    }

    public IBox getBox() {
        return box;
    }

    public void setBox(IBox box) {
        this.box = box;
        if (chainLengthManager != null) {
            // allow old agentManager to de-register itself as a BoxListener
            chainLengthManager.dispose();
        }
        chainLengthManager = new AtomLeafAgentManager(this,box);
    }

    public boolean accept(IAtom a) {
        return ((LengthAgent)chainLengthManager.getAgent((IAtomLeaf)a)).chainLength == maxChainLength;
    }

    private static final long serialVersionUID = 1L;
    protected IBox box;
    protected AtomLeafAgentManager chainLengthManager;
    protected AtomLeafAgentManager agentManager;
    protected int maxChainLength = 0;

    public static class LengthAgent {
        public int chainLength;
    }
}
