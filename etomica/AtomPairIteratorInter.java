package etomica; 

/**
 * Pair iterator for case in which the atom iterators apply to 
 * different sets of atoms.
 *
 * @author David Kofke
 */
public class AtomPairIteratorInter extends AtomPairIterator {
    
    private AtomIteratorSinglet singlet = new AtomIteratorSinglet();
    private AtomIterator ai1Iter, ai1Nbr, ai2Iter, ai2Nbr;
    
    /**
     * @param iter1    iterator of atoms in first group
     * @param iter2    iterator of atoms in second group
     * @param iter1Nbr neighbor-type iterator of atoms in first group
     * @param iter2Nbr neighbor-type iterator of atoms in second group
     */
    public AtomPairIteratorInter(Phase p, AtomIterator iter1, AtomIterator iter2,
                                    AtomIterator iter1Nbr, AtomIterator iter2Nbr) {
        super(p, iter1, iter2Nbr);
        ai1Iter = iter1;
        ai2Iter = iter2;
        ai1Nbr = iter1Nbr;
        ai2Nbr = iter2Nbr;
    }
    
    //to be completed
    public void reset(IteratorDirective id) {}
        
    public void reset() {
        ai1.reset();
        setFirst();
    }
        
    public void reset(Atom atom) {
        if(ai1Iter.contains(atom)) {
            ai1 = ai1Iter;
            ai2 = ai2Nbr;
        }
        else {
            ai1 = ai2Iter;
            ai2 = ai1Nbr;
        }
        ai1.reset(atom);
        setFirst();
    }
        
    public void reset(Atom atom1, Atom atom2) {
        reset(atom1);
        ai1.reset(atom1, atom2);
        setFirst();
    }
}//end of Inter    
        
