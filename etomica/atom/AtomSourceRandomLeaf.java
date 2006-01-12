package etomica.atom;

import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.util.Debug;

/**
 * AtomSource that returns a completely random leaf atom.
 */
public class AtomSourceRandomLeaf implements AtomSource, java.io.Serializable {

    public void setPhase(Phase p) {
        list = p.getSpeciesMaster().leafList;
    }
    
    /**
     * returns a random atom from the phase's leaf atom list
     */
    public Atom getAtom() {
        if (Debug.ON && list== null) throw new IllegalStateException("must set the phase before calling getAtom");
        return list.get(Simulation.random.nextInt(list.size()));
    }
    
    protected AtomArrayList list = null;
}
