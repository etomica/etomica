package etomica.action;

import etomica.api.IAtomLeaf;

/**
 * Interface for a class that can perform an action on an atom.
 */
public interface AtomAction {

    public void actionPerformed(IAtomLeaf atom);
}
