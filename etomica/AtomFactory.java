package etomica;

/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 */
public abstract class AtomFactory implements java.io.Serializable {
    
    public abstract Atom makeAtom(AtomGroup parent, int index);
    
    /**
     * Indicator of whether this factory makes atom groups, or if it 
     * makes atoms that are not groups of other atoms (instance of AtomGroup).
     */
    public abstract boolean producesAtomGroups();
    
    
    public abstract boolean vetoAddition(Atom a); //be sure to check that a is non-null
    
}