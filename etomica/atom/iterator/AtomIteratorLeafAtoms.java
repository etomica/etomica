package etomica.atom.iterator;

import etomica.phase.Phase;

/**
 * Iterator that will loop over all leaf atoms in a phase. Can be configured to
 * iterate all leaf atoms, or only those of a particular species.
 */
public final class AtomIteratorLeafAtoms extends AtomIteratorAdapter implements 
        AtomIteratorPhaseDependent {

    /**
     * Creates iterator with no phase specified. Iteration will return no atoms
     * until a call to setPhase is performed.
     */
    public AtomIteratorLeafAtoms() {
        super(new AtomIteratorArrayListSimple());
    }

    /**
     * Creates iterator conditioned to give all leaf atoms of the specified
     * phase. Call to reset() is required before beginning iteration.
     */
    public AtomIteratorLeafAtoms(Phase phase) {
        this();
        setPhase(phase);
    }

    /**
     * Configures iterator to form its iterates from the leaf atoms of the given
     * phase. If a species was previously (or subsequently) set, iterates will
     * be the leaf atoms of under the species in the specified phase.
     * @throws a NullPointerException if the Phase is null
     */
    public void setPhase(Phase phase) {
        ((AtomIteratorArrayListSimple)iterator).setList(phase.getSpeciesMaster().getLeafList());
    }

    private static final long serialVersionUID = 1L;
}
