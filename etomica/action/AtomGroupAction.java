package etomica.action;

import java.io.Serializable;

import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IMolecule;

/**
 * Wraps an AtomAction, and performs the wrapped action on the atom
 * only if it is a leaf atom; if given an atom group (as indicated
 * by the atom's node), performs action instead on all the atom's
 * child atoms.  This process continues recursively until the leaf atoms are
 * encountered.
 *
 * @author David Kofke
 */
public class AtomGroupAction implements AtomAction, Serializable {

    /**
     * Constructor takes wrapped action, which is final.
     */
    public AtomGroupAction(AtomAction action) {
        this.action = action;
    }
    /* (non-Javadoc)
     * @see etomica.action.AtomAction#actionPerformed(etomica.Atom)
     */
    public void actionPerformed(IAtom atom) {
        if(atom instanceof IMolecule) {
            AtomSet atomList = ((IMolecule)atom).getChildList();
            int size = atomList.getAtomCount();
            for(int i=0; i<size; i++) {
                this.actionPerformed(atomList.getAtom(i));
            }
        }
        else {
            action.actionPerformed(atom);
        }
    }

    /**
     * @return Returns the wrapped action.
     */
    public AtomAction getAction() {
        return action;
    }

    private static final long serialVersionUID = 1L;
    private final AtomAction action;
}
