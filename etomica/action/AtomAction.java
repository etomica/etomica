package etomica.action;

import etomica.Atom;

/**
 * Interface for a class that can perform an action on an atom.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface AtomAction extends AtomsetAction {

    public void actionPerformed(Atom atom);
    public void setAtom(Atom atom);
    public Atom getAtom();
    
    public static final AtomAction NULL = new AtomActionAdapter() {
        public final void actionPerformed(Atom a) {}
    };

}
