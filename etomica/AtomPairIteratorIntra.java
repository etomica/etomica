package etomica; 

/**
 * Pair iterator for case in which the atom iterators apply to the 
 * same set of atoms.
 *
 * @author David Kofke
 */
public class AtomPairIteratorIntra extends AtomPairIterator {
    
    private AtomIteratorSinglet singlet = new AtomIteratorSinglet();
    private AtomIterator ai1Up, ai2Up, ai1Dn, ai2Dn;
    
    public AtomPairIteratorIntra(Phase p, AtomIterator iterUp, AtomIterator iterUpNbr,
                                          AtomIterator iterDn, AtomIterator iterDnNbr) {
        super(p, iterUp, iterUpNbr);
        ai1Up = iterUp;
        ai2Up = iterUpNbr;
        ai1Dn = iterDn;
        ai2Dn = iterDnNbr;
    }
        
    public void reset(IteratorDirective id) {
        if(id.direction() == IteratorDirective.UP) {
            ai1 = ai1Up;
            ai2 = ai2Up;
        }
        else {
            ai1 = ai1Dn;
            ai2 = ai2Dn;
        }
        switch(id.atomCount()) {
            case 0:  ai1.reset(); 
                     break;
            case 1:  ai1 = singlet;
                     ai1.reset(id.atom1());
                     break;
            case 2:  ai1.reset(id.atom1(), id.atom2()); 
                     break;
            default: hasNext = false; 
                     return;
        }//end switch
        setFirst();
    }
}//end of Intra    
        
