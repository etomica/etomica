package etomica.atom;


/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 *
 * @author David Kofke and Andrew Schultz
 */
 
public abstract class AtomFactory implements java.io.Serializable {
    
    protected final AtomType atomType;
    protected boolean isMutable;
    
    public AtomFactory(AtomType atomType) {
        this.atomType = atomType;
        isMutable = true;
    }
    
    /**
     * Builds and returns the atom/atomgroup made by this factory.
     * Implementation of this method in the subclass defines the 
     * product of this factory.
     */
    public abstract IAtom makeAtom();

    /**
     * Returns the number of number of atoms used to form the Atom
     * returned by makeAtom.  This includes the Atom itself, its children,
     * their children, etc., down to the leaf atoms.
     */
    public abstract int getNumTreeAtoms();
    
    /**
     * Returns the number of child atoms held by the Atom returned by
     * makeAtom.  This will be zero if makeAtom returns a leaf atom.
     */
    public abstract int getNumChildAtoms();
    
    /**
     * Returns the number of leaf atoms descended from the Atom returned 
     * by makeAtom.  This will be 1 if makeAtom returns a leaf atom.
     */
    public abstract int getNumLeafAtoms();
    
    /**
     * Returns the atomType instance given to all atoms made by this factory.
     */
    public AtomType getType() {
        return atomType;
    }
    
    public boolean isMutable() {
        return isMutable;
    }
}
