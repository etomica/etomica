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
public class AtomPair implements AtomSet {


    /* (non-Javadoc)
     * @see etomica.AtomSet#getAtom(int)
     */
    public final Atom getAtom(int i) {
        if(i == 0) return atom0;
        if(i == 1) return atom1;
        throw new IllegalArgumentException();
    }

    /* (non-Javadoc)
     * @see etomica.AtomSet#count()
     */
    public final int count() {
        return 2;
    }
    
    public void copyTo(AtomPair pair) {
        pair.atom0 = atom0;
        pair.atom1 = atom1;
    }
    
    public boolean equals(AtomSet pair) {
        return pair instanceof AtomPair
                && ((AtomPair)pair).atom0 == atom0
                && ((AtomPair)pair).atom1 == atom1;
    }
    
    public String toString() {
        return atom0+" and "+atom1;
    }
    
    public Atom atom0, atom1;

}
