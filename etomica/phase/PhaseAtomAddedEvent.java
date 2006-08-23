package etomica.phase;

import etomica.atom.Atom;


/**
 * Event that conveys that an Atom has been added to a Phase.
 */
public class PhaseAtomAddedEvent extends PhaseAtomEvent {

    public PhaseAtomAddedEvent(Phase phase, Atom atom) {
        super(phase, atom);
    }

    private static final long serialVersionUID = 1L;
}
