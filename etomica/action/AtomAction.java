package etomica.action;

import etomica.Action;
import etomica.Atom;

/**
 * @author David Kofke
 *
 * Interface for a class that can perform an action on an atom.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public interface AtomAction extends Action {

    public void actionPerformed(Atom atom);
    public void setAtom(Atom atom);
    public Atom getAtom();
    
    public static final AtomActionAdapter NULL = new Null();

    /**
     * Action that does nothing.
     */
    public static class Null extends AtomActionAdapter {
        public final void actionPerformed(Atom a) {}
    }

}
