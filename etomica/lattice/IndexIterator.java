package etomica.lattice;


/**
 * Loops through the set of indexes appropriate to a lattice of given size.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jan 3, 2005 by kofke
 */
public interface IndexIterator {

    public void reset();
    
    public boolean hasNext();
    
    public int[] next();
}
