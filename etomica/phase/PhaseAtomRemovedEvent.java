package etomica.phase;

import etomica.atom.IAtom;


/**
 * Event that conveys that an Atom has been removed from a Phase.
 */
public class PhaseAtomRemovedEvent extends PhaseAtomEvent {

    public PhaseAtomRemovedEvent(Phase phase, IAtom atom) {
        super(phase, atom);
    }

    private static final long serialVersionUID = 1L;
}
