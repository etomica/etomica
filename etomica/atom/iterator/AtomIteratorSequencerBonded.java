package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;


/**
 * 
 */

//TODO fix this class
public class AtomIteratorSequencerBonded extends AtomIteratorSequencerList {

    public AtomIteratorSequencerBonded() {
        super();
//        throw new RuntimeException("this class needs to be updated to handle direction correctly");
    }
    public boolean contains(AtomSet atom) {
        reset();
        return listIterator.peek() == atom; 
    }

    public int size() {
        reset();
        return listIterator.hasNext() ? 1 : 0;
    }

    public Atom nextAtom() {
        Atom next = listIterator.nextAtom();
        listIterator.unset();
        return next;
    }

    public AtomSet next() {
        return nextAtom();
    }
    

}
