package etomica.phase;

import etomica.atom.IAtom;


/**
 * Event that conveys that an Atom's global index in a phase has changed.
 */
public class PhaseAtomIndexChangedEvent extends PhaseAtomEvent {

    public PhaseAtomIndexChangedEvent(Phase phase, IAtom atom, int oldIndex) {
        super(phase, atom);
        this.oldIndex = oldIndex;
    }
    
    public int getOldIndex() {
        return oldIndex;
    }
    
    private final int oldIndex;
    private static final long serialVersionUID = 1L;
}
