package etomica;

/**
 * Forms pairs from 1 atom of one group with All the atoms from another group.
 *
 * @author David Kofke
 */
public final class ApiIntergroup1A implements AtomPairIterator {
    
    public ApiIntergroup1A(Simulation sim) {
        pair = new AtomPair(sim.space);
        iterator = new AtomIteratorListSimple();
    }

    public void setBasis(Atom a1, Atom a2) {
        if(a1 == a2 || a1 == null || a2 == null) throw new IllegalArgumentException("Improper basis given to ApiIntergroup1A");
        basis1 = a1;
        basis2 = a2;
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator, based
     * on the most recent specification of the atom given to reset.
     */
    public int size() {return iterator.size();}       
    
    public boolean hasNext() {return hasNext;}
    
    /**
     * Uses atom in iterator directive to specify the single atom to be given
     * in all iterated pairs.  This atom will be the child of the
     * current basis that has the specified atom as a descendant.
     */
    public void reset(IteratorDirective id) {
        reset(id.atom1());
    }
    
    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        iterator.reset();
        hasNext = (atom1 != null && iterator.hasNext());
        pair.atom1 = atom1;
    }
        
    /**
     * Uses given atom to specify the single atom to be present
     * in all iterated pairs.  This single atom will be the child of the
     * current basis that has the given atom as a descendant.
     */
    public void reset(Atom atom) {
        atom1 = atom;
        if(atom1 == null) {hasNext = false; return;}
        
        //find which basis this atom descends from, then set iterator to other basis
        //move up atom tree until one basis is encountered
        while(atom1 != null) {
            Atom parent = atom1.node.parentGroup();
            if(parent == basis1) {
                iterator.setBasis(basis2);
                break;
            } else if(parent == basis2) {
                iterator.setBasis(basis1);
                break;
            }
            atom1 = parent;
        }
        reset();
    }
        
    public AtomPair next() {
        pair.atom2 = iterator.next();
        pair.reset();
        hasNext = iterator.hasNext();
        return pair;
    }

     //needs to change for neighbor iteration
    public void allPairs(AtomPairAction act) {
        throw new RuntimeException("ApiIntergroup1A.allPairs not yet implemented");
   /*     Atom basis = iterator.getBasis();
        if(basis == null || atom1 == null) return;
        Atom last = basis.node.lastChildAtom();
        for(Atom atom = basis.node.firstChildAtom(); atom != null; atom=atom.nextAtom()) {
            pair.atom2 = atom;
            pair.reset();
            act.action(pair);
            if(atom == last) break;
        }*/
    }
    
    private final AtomIteratorListSimple iterator;
    private Atom basis1, basis2;
    private Atom atom1;
    private boolean hasNext;
    private final IteratorDirective localDirective = new IteratorDirective();
    
    private final AtomPair pair;
    
}  //end of class AtomPairIterator
    
