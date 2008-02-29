package etomica.action;

import java.io.Serializable;

import etomica.api.IVector;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.space.Space;

/**
 * 
 * Moves (translates) an atom by a specified vector amount.
 * To move all atoms in a molecule (or atom group), wrap an
 * instance of this class in an AtomGroupAction.
 * 
 * @author David Kofke
 */
public class AtomActionTranslateBy implements AtomAction, Serializable {
    
    private static final long serialVersionUID = 1L;
    private final IVector translationVector;
    
    public AtomActionTranslateBy(Space space) {
        translationVector = space.makeVector();
    }
    
    public void actionPerformed(IAtom atom) {
        ((IAtomPositioned)atom).getPosition().PE(translationVector);
    }
       
    /**
     * Returns the translation vector, the distance that the
     * atom will be moved by this action. Returns the vector used by this
     * instance, not a copy, so any manipulation of the returned vector will
     * affect the action of this instance.
     */
    public IVector getTranslationVector() {
        return translationVector;
    }
    /**
     * @param destination The translation vector to set.  A local copy
     * is made of the given vector.
     */
    public void setTranslationVector(IVector translationVector) {
        this.translationVector.E(translationVector);
    }
}