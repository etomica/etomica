package etomica;

/**
 * Interface that defines the default beginning and end points of 
 * iteration for an atom iterator.
 *
 * @author David Kofke
 */
 
 public interface AtomIteratorBasis {
    
    public Atom firstAtom();
    
    public Atom lastAtom();
    
    public boolean contains(Atom a);
 }