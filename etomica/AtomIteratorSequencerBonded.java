package etomica;

import etomica.action.AtomsetAction;

/**
 * 
 */
public class AtomIteratorSequencerBonded extends AtomIteratorSequencerList {

    public void allAtoms(AtomsetAction action) {
        reset();
        if (listIterator.hasNext()) {
            action.actionPerformed(listIterator.next());
        }
    }

    public boolean contains(Atom[] atom) {
        return listIterator.peek()[0] == atom[0]; 
    }

    public int size() {
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
