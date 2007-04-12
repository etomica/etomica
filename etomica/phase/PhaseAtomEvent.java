package etomica.phase;

import etomica.atom.IAtom;

/**
 * Event that conveys some happening with respect to an Atom in a Phase.
 */
public class PhaseAtomEvent extends PhaseEvent {
    
    public PhaseAtomEvent(Phase phase, IAtom atom) {
        super(phase);
        this.atom = atom;
    }

    public IAtom getAtom() {
        return atom;
    }
    
    private final IAtom atom;
    private static final long serialVersionUID = 1L;
}
