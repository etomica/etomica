package etomica.atom;

import etomica.box.Box;

public class AtomSourceRandomMoleculeSeq extends AtomSourceRandomMolecule {

    public void setBox(Box p) {
        super.setBox(p);
        atomManager = p.getSpeciesMaster();
        agentList = atomManager.getAgentList();
        reset();
    }
    
    public IAtom getAtom() {
        int moleculeCount = atomManager.moleculeCount();
        // this is probably innapropriate for atom removal
        if (prevIndex == -1 || prevIndex > moleculeCount-1) {
            // no suitable previous atom to step forward from the beginning
            reset();
            // no atoms in the box
            if (prevIndex == -1) return null;
        }
        int lookAhead = random.nextInt(maxLookAhead+1);
        
        for ( ; agentIndex<agentList.getAtomCount(); agentIndex++) {
            // advance through the species if needed
            moleculeList = ((IAtomGroup)agentList.getAtom(agentIndex)).getChildList();
            int count = moleculeList.getAtomCount();
            if (prevIndex+lookAhead < count) {
                prevIndex += lookAhead;
                lookAhead = -1;
            }
            prevIndex -= count;
        }
        if (lookAhead > -1) {
            // we ran out of species, so start over with the first species
            for (agentIndex=0 ; agentIndex<agentList.getAtomCount(); agentIndex++) {
                moleculeList = ((IAtomGroup)agentList.getAtom(agentIndex)).getChildList();
                int count = moleculeList.getAtomCount();
                if (prevIndex+lookAhead < count) {
                    prevIndex += lookAhead;
                    lookAhead = -1;
                }
                prevIndex -= count;
            }
        }
        return moleculeList.getAtom(prevIndex);
    }

    /**
     * Reset the atom used to step from to a random molecule
     */
    public void reset() {
        int size = atomManager.moleculeCount();
        if (size == 0) {
            prevIndex = -1;
            return;
        }
        prevIndex = random.nextInt(size);
        
        for (agentIndex=0; agentIndex<agentList.getAtomCount(); agentIndex++) {
            moleculeList = ((IAtomGroup)agentList.getAtom(agentIndex)).getChildList();
            int count = moleculeList.getAtomCount();
            if (prevIndex < count) {
                break;
            }
            prevIndex -= count;
        }
    }

    /**
     * Returns the maximum number of molecules the source will advance for
     * each call to getAtom
     */
    public int getMaxLookAhead() {
        return maxLookAhead;
    }
    /**
     * Sets the maximum number of molecules the source will advance for
     * each call to getAtom
     */
    public void setMaxLookAhead(int maxLookAhead) {
        this.maxLookAhead = maxLookAhead;
    }
    
    private static final long serialVersionUID = 2L;
    protected int maxLookAhead = 10;
    protected AtomManager atomManager;
    protected AtomSet moleculeList, agentList;
    protected int prevIndex, agentIndex;
}