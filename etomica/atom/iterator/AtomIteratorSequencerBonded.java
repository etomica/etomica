package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;


/**
 * 
 */

//TODO fix this class
public class AtomIteratorSequencerBonded extends AtomIteratorSequenceDirectable {

    public AtomIteratorSequencerBonded() {
        super();
//        throw new RuntimeException("this class needs to be updated to handle direction correctly");
    }
    public boolean contains(AtomSet atom) {
        reset();
        return iterator.peek() == atom; 
    }

    public int size() {
        reset();
        return iterator.hasNext() ? 1 : 0;
    }

    public Atom nextAtom() {
        Atom next = iterator.nextAtom();
        if (doGoDown)
            resetDown();
        else
            iterator.unset();
        return next;
    }

}
