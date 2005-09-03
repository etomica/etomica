package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomSet;
import etomica.atom.iterator.IteratorDirective.Direction;


/**
 * Returns one or both of the atoms adjacent to a specified atom in
 * its sequence list (i.e., the AtomList containing its seq linker).
 * If adjacent linker has no atom (is a tab or header) no corresponding
 * atom is given as an iterate; thus iterator may give 0, 1, or 2 iterates
 * depending on presence of adjacent atoms and specification of iteration
 * direction.
 */

public class AtomIteratorSequenceAdjacent implements AtomIteratorAtomDependent, java.io.Serializable {

    /**
     * Constructor gives iterator not ready for iteration.  Must
     * set an atom and reset before using.  Default direction is
     * null, meaning that both adjacent atoms (if there are two)
     * will be given by iterator.
     */
    public AtomIteratorSequenceAdjacent(IteratorDirective.Direction direction) {
        super();
        this.direction = direction;
        setAtom(null);
        unset();
    }

    /**
     * Returns true if the given AtomSet has count == 1 and 
     * its atom is one of the iterates for the current condition
     * of the iterator (independent of hasNext status).
     */
    public boolean contains(AtomSet atom) {
        if(atom == null || atom.count() != 1) return false;
        Atom testAtom = atom.getAtom(0);
        if(testAtom == null) return false;
        if(direction == IteratorDirective.UP && (atomSeq.next.atom == testAtom)) {
            return true;
        }
        if(direction == IteratorDirective.DOWN && (atomSeq.previous.atom == testAtom)) {
            return true;
        }
        return false;
    }

    /**
     * Returns the number of iterates that iterator would give
     * if reset and iterated in its current condition.  Does not
     * depend on or affect iteration state.
     */
    public int size() {
        if(direction == IteratorDirective.UP && (atomSeq.next.atom != null)) {
            return 1;
        }
        if(direction == IteratorDirective.DOWN && (atomSeq.previous.atom != null)) {
            return 1;
        }
        return 0;
    }

    /**
     * Performs action on all iterates for current condition of iterator.
     */
    public void allAtoms(AtomsetAction action) {
        if(direction == IteratorDirective.UP) {
            Atom atom = atomSeq.next.atom;
            if(atom != null) action.actionPerformed(atom);
        }
        else {
            Atom atom = atomSeq.previous.atom;
            if(atom != null) action.actionPerformed(atom);
        }
    }

    /**
     * Returns true if another iterate is forthcoming, false otherwise.
     */
    public boolean hasNext() {
        return hasNext;
    }

    /**
     * Returns 1, indicating that this is an atom AtomSet iterator.
     */
    public int nBody() {
        return 1;
    }

    /**
     * Same as nextAtom.
     */
    public AtomSet next() {
        return nextAtom();
    }

    /**
     * Returns the next iterator, or null if hasNext is false.
     */
    public Atom nextAtom() {
        if (!hasNext) {
            return null;
        }
        hasNext = false;
        if(direction == IteratorDirective.UP) {
            return atomSeq.next.atom;
        }
        return atomSeq.previous.atom;
    }

    /**
     * Returns the next iterate without advancing the iterator.
     */
    public AtomSet peek() {
        if (!hasNext) {
            return null;
        }
        if(direction == IteratorDirective.UP) {
            return atomSeq.next.atom;
        }
        return atomSeq.previous.atom;
    }

    /**
     * Readies the iterator to begin iteration.
     */
    public void reset() {
        hasNext = ((direction != IteratorDirective.DOWN) && (atomSeq.next.atom != null)) ||
                  ((direction != IteratorDirective.UP) && (atomSeq.previous.atom != null));
    }

    /**
     * Puts iterator in a state where hasNext is false.
     */
    public void unset() {
        hasNext = false;
    }   
    
    /**
     * Sets the central atom.  Iterates will be those atoms (if not null) at
     * atom.seq.next.atom and atom.seq.previous.atom, as indicated
     * by the direction field.
     */
    public void setAtom(Atom atom) {
        atomSeq = (atom != null) ? atom.seq : emptyList.header;
    }
    
    private AtomLinker atomSeq;
    private final Direction direction;
    private final AtomList emptyList = new AtomList();
    private boolean hasNext;
}
