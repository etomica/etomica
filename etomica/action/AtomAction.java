package etomica.action;

import etomica.atom.Atom;

/**
 * Interface for a class that can perform an action on an atom.
 */
public interface AtomAction extends AtomsetAction {

    public void actionPerformed(Atom atom);
    public void setAtom(Atom atom);
    public Atom getAtom();
    
}
