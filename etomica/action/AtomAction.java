package etomica.action;

import etomica.atom.IAtom;

/**
 * Interface for a class that can perform an action on an atom.
 */
public interface AtomAction {

    public void actionPerformed(IAtom atom);
}
