package etomica.phase;

import etomica.atom.Atom;

/**
 * Event that conveys some happening with respect to an Atom in a Phase.
 */
public class PhaseAtomEvent extends PhaseEvent {
    
    public PhaseAtomEvent(Phase phase, Atom atom) {
        super(phase);
        this.atom = atom;
    }

    public Atom getAtom() {
        return atom;
    }
    
    private final Atom atom;
    private static final long serialVersionUID = 1L;
}
