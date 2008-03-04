package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomType;



/**
 * Filter that accepts only those atoms have a specific type.  
 * Types must be the same instance to match; different instances
 * of the same type are not accepted.
 *
 * @author David Kofke
 *
 */

public class AtomFilterTypeInstance implements AtomFilter, java.io.Serializable {

    /**
     * @param acceptedType type instance that must equal atom type for atom
     * to be accepted 
     */
    public AtomFilterTypeInstance(IAtomType acceptedType) {
        this.acceptedType = acceptedType;
    }

    /**
     * Returns true if the type of the given atom is the same 
     * instance as the type given at construction.  Returns false
     * also if atom is null.
     */
    public boolean accept(IAtom atom) {
        return (atom != null) && (atom.getType() == acceptedType);
    }

    /**
     * @return Returns the acceptedType.
     */
    public IAtomType getAcceptedType() {
        return acceptedType;
    }
    
    private static final long serialVersionUID = 1L;
    private final IAtomType acceptedType;
}
