package etomica.action;
import etomica.Atom;
import etomica.Space;
import etomica.space.Vector;

/**
 * 
 * Moves (translates) an atom by a specified vector amount.
 * To move all atoms in a molecule (or atom group), wrap an
 * instance of this class in an AtomGroupAction.
 * 
 * @author David Kofke
 */
public class AtomActionTranslateBy extends AtomActionAdapter {
    
    private final Vector translationVector;
    
    public AtomActionTranslateBy(Space space) {
        translationVector = space.makeVector();
    }
    
    public void actionPerformed(Atom atom) {
        atom.coord.position().PE(translationVector);
    }
       
    /**
     * Returns the translation vector, the distance that the
     * atom will be moved by this action. Returns the vector used by this
     * instance, not a copy, so any manipulation of the returned vector will
     * affect the action of this instance.
     */
    public Vector getTranslationVector() {
        return translationVector;
    }
    /**
     * @param destination The translation vector to set.  A local copy
     * is made of the given vector.
     */
    public void setTranslationVector(Vector translationVector) {
        this.translationVector.E(translationVector);
    }
}