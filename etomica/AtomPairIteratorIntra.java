package etomica; 

/**
 * Pair iterator for case in which the atom iterators apply to the 
 * same set of atoms.
 *
 * @author David Kofke
 */
public class AtomPairIteratorIntra extends AtomPairIterator {
    
    private AtomIterator.Singlet singlet = new AtomIterator.Singlet();
    private AtomIterator ai1Save;
    
    public AtomPairIteratorIntra(Phase p, AtomIterator iter, AtomIterator iterNbr) {
        super(p, iter, iterNbr);
        ai1Save = ai1;
    }
        
    public void reset() {
        ai1 = ai1Save;
        ai1.reset();
        setFirst();
    }
        
    public void reset(Atom atom) {
        ai1 = singlet;
        ai1.reset(atom);
        setFirst();
    }
        
    public void reset(Atom atom1, Atom atom2) {
        ai1 = ai1Save;
        ai1.reset(atom1, atom2);
        setFirst();
    }
}//end of Intra    
        
