package etomica;

/**
 * Class responsible for building new instances of the atoms (or atom groups)
 * that are collected in a given AtomGroup.
 */
public abstract class AtomFactory implements java.io.Serializable {
    
    public abstract Atom makeAtom(AtomGroup parent, int index);
    
}