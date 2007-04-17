package etomica.atom.iterator;

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
     * @throws a NullPointerException if the Phase is null
     */
    public void setPhase(Phase phase) {
        ((ApiIntraArrayList)iterator).setList(phase.getSpeciesMaster().getLeafList());
    }

    private static final long serialVersionUID = 1L;
}
