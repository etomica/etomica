package etomica;

/**
 * Interface for class that iterates over triples of atoms.
 *
 * @author David Kofke
 */
public interface Atom3Iterator extends AtomSetIterator, java.io.Serializable {
    
    public boolean hasNext();
    public void reset(IteratorDirective id);

    /**
     * Resets the iterator, so that it is ready to go through all of its triplets.
     */
    public void reset();
    
    public Atom3 next();
    
    public void setBasis(Atom a1, Atom a2, Atom a3);

}  //end of interface Atom3Iterator
    
