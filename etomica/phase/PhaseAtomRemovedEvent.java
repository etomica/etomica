package etomica.phase;

import etomica.atom.Atom;


/**
 * Event that conveys that an Atom has been removed from a Phase.
 */
public class PhaseAtomRemovedEvent extends PhaseAtomEvent {

    public PhaseAtomRemovedEvent(Phase phase, Atom atom) {
        super(phase, atom);
    }

    private static final long serialVersionUID = 1L;
}
