package etomica;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Feb 18, 2005 by kofke
 */
public interface AtomSet {

    /**
     * Returns the i-th atom, with numbering beginning from 0. 
     * If i is greater than count-1, throws an IllegalArgumentException.
     */
    public Atom getAtom(int i);
    
    public int count();
    
    public boolean equals(AtomSet atoms);
    
}
