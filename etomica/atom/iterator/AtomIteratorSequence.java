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
 * with a specified linker. This linker may be specified directly, or as one
 * associated with a specified atom via an AtomToLinker instance (given at
 * construction - default uses the atom's sequencer as the linker). Direct,
 * unterminating iteration of linkers can be performed using nextLinker method. 
 */

/*
 * Created on June 5, 2005 by kofke
 */
public class AtomIteratorSequence implements AtomIteratorAtomDependent, java.io.Serializable {

    /**
     * Constructs a new class to iterate in the specified direction. Default
     * uses the atom's sequencer as the linker to begin atom-specified
     * iteration.
     */
    public AtomIteratorSequence(IteratorDirective.Direction direction) {
        this(direction, DEFAULT);
    }

    /**
     * Constructs new class to iterate in the specified direction. Must invoke
     * setAtom or setFirst and reset before beginning iteration. AtomToLinker
     * instance provides a rule that identifies the linker for beginning
     * iteration, given an atom.
     * 
     * @throws IllegalArgumentException
     *             if direction is null
     */
    public AtomIteratorSequence(IteratorDirective.Direction direction,
            AtomToLinker atomToLinker) {
        if (direction == null)
            throw new IllegalArgumentException(
                    "Must specify direction to constructor of AtomLinkerIterator");
        upListNow = (direction == IteratorDirective.UP);
        this.atomToLinker = atomToLinker;
        setFirst(emptyList.header);
    }

    /**
     * Indicates if iterator has another atom iterate. Has no connection to
     * nextLinker method.
     */
    public boolean hasNext() {
        return next.atom != null;
    }

    /**
     * Same as nextAtom, for AtomsetIterator interface.
     */
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
        next = upListNow ? next.next : next.previous;
        return upListNow ? next.previous : next.next;
    }
    
    /**
     * Returns the next linker without advancing the iterator.
     */
    public AtomLinker peekLinker() {
        return next;
    }

    /**
     * Returns the next atom without advancing the iterator.
     */
    public AtomSet peek() {
        return next.atom;
    }

    /**
     * Puts iterator in a state in which hasNext is false.
     */
    public void unset() {
        next = emptyList.header;
    }

    /**
     * Initializes iterator for iteration.
     */
    public void reset() {
        next = first;
    }

    /**
     * Performs action on all atoms for current condition (as given via setFirst
     * or setAtom) of iterator. Unaffected by and does not affect state (hasNext
     * status) of iterator.
     */
    public void allAtoms(AtomsetAction action) {
        for (AtomLinker link = first; link.atom != null; link = upListNow ? link.next
                : link.previous) {
            action.actionPerformed(link.atom);
        }
    }

    /**
     * Returns true if given atom set is among the iterates, if iterated after a
     * call to reset().
     */
    public boolean contains(AtomSet atom) {
        if (atom == null || atom.count() != 1)
            return false;
        detector.setAtoms(atom);
        detector.reset();
        allAtoms(detector);
        return detector.detectedAtom();
    }

    /**
     * Returns the number of iterates given by this iterator, if iterated after
     * a call to reset().
     */
    public int size() {
        counter.reset();
        allAtoms(counter);
        return counter.callCount();
    }

    /**
     * Returns 1, indicating that this is an Atom iterator.
     */
    public int nBody() {
        return 1;
    }

    /**
     * Sets the first atom for iteration. Iteration proceeds from this atom up
     * and/or down the list, depending on how iterator was configured at
     * construction, and ends when a null-atom linker is encountered.
     */
    public void setAtom(Atom atom) {
        setFirst(atomToLinker.getLinker(atom));
    }

    /**
     * Sets the linker where iteration begins. The given linker's atom will be
     * the first iterate given. If linker has a null atom, no iterates will be
     * given. Otherwise, iteration proceeds from this linker up and/or down the
     * list, depending on how iterator was configured at construction, and
     * ends when a null-atom linker is encountered.
     *  
     */
    public void setFirst(AtomLinker newFirst) {
        first = (newFirst != null) ? newFirst : emptyList.header;
        unset();
    }

    /**
     * @return the linker given at the most recent call to
     *         setFirst, or an empty-list list header if a first linker
     *         specified, or unset has been invoked.
     */
    public AtomLinker getFirst() {
        return first;
    }

    /**
     * Returns the direction of iteration for this iterator as configured
     * at construction.
     */
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
     * Interface for class that determines an atom linker given an atom to begin
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
