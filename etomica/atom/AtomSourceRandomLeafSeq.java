package etomica.atom;

import etomica.phase.Phase;
import etomica.simulation.Simulation;

/**
 * Returns a semi-random leaf atom from a phase.  Each atom is 0 to N atoms 
 * from the previously returned atom.
 */
public class AtomSourceRandomLeafSeq extends AtomSourceRandomLeaf {

    public void setPhase(Phase p) {
        list = p.getSpeciesMaster().atomList;
        reset();
    }
    
    public Atom getAtom() {
        if (prevLinker == null || !prevLinker.atom.node.inSimulation()) {
            // no suitable previous atom to step forward from 
            reset();
            // no atoms in the phase
            if (prevLinker == null) return null;
        }
        int lookAhead = Simulation.random.nextInt(maxLookAhead+1);
        for (int i=0; i<lookAhead; i++) {
            do {
                prevLinker = prevLinker.next;
            } while (prevLinker.atom == null);
        }
        return prevLinker.atom;
    }

    /**
     * Reset the atom used to step from to a random leaf atom
     */
    public void reset() {
        int size = list.size();
        if (size == 0) {
            prevLinker = null;
        }
        else {
            prevLinker = list.entry(Simulation.random.nextInt(size));
        }
    }
        
    /**
     * Returns the maximum number of atoms the AtomSource will step
     * ahead from the previous one.
     */
    public int getMaxLookAhead() {
        return maxLookAhead;
    }

    /**
     * Sets the maximum number of atoms the AtomSource will step
     * ahead from the previous one.
     */
    public void setMaxLookAhead(int maxLookAhead) {
        this.maxLookAhead = maxLookAhead;
    }
    
    protected int maxLookAhead = 10;
    protected AtomLinker prevLinker;
}