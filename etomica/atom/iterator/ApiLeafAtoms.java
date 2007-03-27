/*
 * History
 * Created on Sep 15, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.atom.AtomArrayList;
import etomica.phase.Phase;

/**
 * Iterator that returns all pairs that can be formed from all leaf atoms of a
 * given phase. Wraps an ApiIntraList instance.
 */
public class ApiLeafAtoms extends AtomPairIteratorAdapter implements
        AtomsetIteratorPhaseDependent {

    /**
     * Creates new pair iterator that requires reset() before beginning
     * iteration.
     */
    public ApiLeafAtoms() {
        super(new ApiIntraArrayList());
    }

    /**
     * Conditions iterator to return all leaf-atom pairs from the given phase.
     * If phase is null, no iterates are given.
     */
    public void setPhase(Phase phase) {
        if (phase == null) {
            emptyList.clear();
            ((ApiIntraArrayList)iterator).setList(emptyList);
        } else {
            ((ApiIntraArrayList)iterator).setList(phase.getSpeciesMaster().getLeafList());
        }
    }

    private static final long serialVersionUID = 1L;
    private final AtomArrayList emptyList = new AtomArrayList();
}
