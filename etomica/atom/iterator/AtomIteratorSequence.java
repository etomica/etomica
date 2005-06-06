package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;
import etomica.IteratorDirective;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;

/**
 * Iterator for looping through a sequence of AtomLink-ed atoms relative to a
 * specified atom. Can be configured at construction to loop upList or downList,
 * and direction cannot be changed afterward. Iteration ends at the first
 * encounter of a linker having a null atom (i.e., at any Tab or a list header).
 * No list specification is provided; rather iteration is performed beginning
 * with a linker associated with the atom via an AtomToLinker instance given at
 * construction. Default version uses the atom's sequencer as the linker.
 */

/*
 * Created on June 5, 2005 by kofke
 */
public class AtomIteratorSequence implements AtomIteratorAtomDependent {

    /**
     * Default instance uses the atom's sequencer as the linker to begin list
     * iteration.
     */
    public AtomIteratorSequence(IteratorDirective.Direction direction) {
        this(direction, DEFAULT);
    }

    /**
     * Constructs new class with hasNext as false. Must invoke setAtom and reset
     * before beginning iteration. Argument is class that identifies the linker
     * for beginning iteration, given an atom.
     */
    //see etomica.nbr.cell.ApiIntraspecies1ACell for non-default use of this
    // constructor
    public AtomIteratorSequence(IteratorDirective.Direction direction,
            AtomToLinker atomToLinker) {
        if (direction == null)
            throw new IllegalArgumentException(
                    "Must specify direction to constructor of AtomLinkerIterator");
        upListNow = (direction == IteratorDirective.UP);
        this.atomToLinker = atomToLinker;
        setFirst(emptyList.header);
    }

    public boolean hasNext() {
        return next.atom != null;
    }

    public final AtomSet next() {
        return nextAtom();
    }

    /**
     * Returns the next atom, and advances iterator if the atom is not null.
     * Repeated calls to this method will not progress beyond any tabs or header
     * encountered in the iteration.
     */
    public Atom nextAtom() {
        Atom atom = next.atom;
        if (atom != null) {
            next = upListNow ? next.next : next.previous;
        }
        return atom;
    }

    /**
     * Returns the next linker and advances iterator, regardless of hasNext
     * status. Repeated calls to this method would loop endlessly through a loop
     * of linkers, returning any tabs, headers, or regular linkers encountered
     * during the iteration.
     */
    public AtomLinker nextLinker() {
        AtomLinker thisLinker = next;
        next = upListNow ? next.next : next.previous;
        return thisLinker;
    }//end of nextLinker

    public AtomSet peek() {
        return next.atom;
    }

    public void unset() {
        next = emptyList.header;
    }

    public void reset() {
        next = first;
    }

    public void allAtoms(AtomsetAction action) {
        for (AtomLinker link = first; link.atom != null; 
                link = upListNow ? link.next : link.previous) {
            action.actionPerformed(link.atom);
        }
    }

    public boolean contains(AtomSet atom) {
        detector.setAtoms(atom);
        detector.reset();
        allAtoms(detector);
        return detector.detectedAtom();
    }

    public int size() {
        counter.reset();
        allAtoms(counter);
        return counter.callCount();
    }

    public int nBody() {
        return 1;
    }

    /**
     * Sets the first atom for iteration. Iteration proceeds from this atom up
     * and/or down the list, as specified by setDirection and setNumToSkip.
     * Atom's sequencer is used to identify its position in the list.
     */
    public void setAtom(Atom atom) {
        setFirst(atomToLinker.getLinker(atom));
    }

    public void setFirst(AtomLinker newFirst) {
        first = (newFirst != null) ? newFirst : emptyList.header;
        unset();
    }

    /**
     * @return the linker of the atom corresponding to the most recent call to
     *         setFirst, or the list header if a first atom has not be
     *         specified.
     */
    public AtomLinker getFirst() {
        return first;
    }

    public IteratorDirective.Direction getDirection() {
        return upListNow ? IteratorDirective.UP : IteratorDirective.DOWN;
    }

    private final boolean upListNow;
    private final AtomToLinker atomToLinker;
    private final AtomsetCount counter = new AtomsetCount();
    private final AtomsetDetect detector = new AtomsetDetect(null);
    private final AtomList emptyList = new AtomList();
    private AtomLinker next;//nonnull
    private AtomLinker first;//nonnull

    /**
     * Interface for class the determines the atom linker given an atom to begin
     * iteration.
     */
    public interface AtomToLinker {

        /**
         * Returns a linker that this instance associates with the given atom.
         * Should return null if the atom is null.
         */
        public AtomLinker getLinker(Atom atom);
    }

    private static final AtomToLinker DEFAULT = new AtomToLinker() {

        public AtomLinker getLinker(Atom atom) {
            return (atom != null) ? atom.seq : null;
        }
    };
}