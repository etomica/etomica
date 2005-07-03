/*
 * History
 * Created on Sep 15, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Phase;
import etomica.atom.AtomList;

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
        super(new ApiIntraList());
        apiList = (ApiIntraList) iterator;
    }

    /**
     * Conditions iterator to return all leaf-atom pairs from the given phase.
     * If phase is null, no iterates are given.
     */
    public void setPhase(Phase phase) {
        if (phase == null) {
            emptyList.clear();
            apiList.setList(emptyList);
        } else {
            apiList.setList(phase.speciesMaster().atomList);
        }
    }

    private final ApiIntraList apiList;
    private final AtomList emptyList = new AtomList();
}
