package etomica.phase;

import etomica.atom.IAtom;


/**
 * Event that conveys that an Atom has been added to a Phase.
 */
public class PhaseAtomAddedEvent extends PhaseAtomEvent {

    public PhaseAtomAddedEvent(Phase phase, IAtom atom) {
        super(phase, atom);
    }

    private static final long serialVersionUID = 1L;
}
