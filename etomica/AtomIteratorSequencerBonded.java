package etomica;


/**
 * 
 */

//TODO fix this class
public class AtomIteratorSequencerBonded extends AtomIteratorSequencerList {

    public AtomIteratorSequencerBonded() {
        super();
        throw new RuntimeException("this class needs to be updated to handle direction correctly");
    }
    public boolean contains(Atom[] atom) {
        reset();
        return listIterator.peek()[0] == atom[0]; 
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

    public Atom[] next() {
        Atom[] next = listIterator.next();
        listIterator.unset();
        return next;
    }
    

}
