package etomica;


/**
 * 
 */
public class AtomIteratorSequencerBonded extends AtomIteratorSequencerList {

    public AtomIteratorSequencerBonded() {
        super();
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
        Atom[] next = iterator.next();
        listIterator.unset();
        return next;
    }
    

}
