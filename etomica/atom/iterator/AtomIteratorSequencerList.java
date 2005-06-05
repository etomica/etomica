package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;
import etomica.IteratorDirective;
import etomica.action.AtomsetAction;
import etomica.atom.AtomLinker;

/**
 * Iterator for looping through the sequence list relative to a specified atom.
 * Can be configured to loop upList, downList, or both direction from specified
 * atom. Also can be configured to skip any number of atoms from the specified
 * atom. No list specification is provided; rather iteration is performed
 * beginning with a linker associated with the atom via a AtomToLinker class
 * given at construction. Default version uses the atom's sequencer as the
 * linker.
 */

/*
 * History Created on Aug 26, 2004 by kofke
 */
public class AtomIteratorSequencerList implements AtomIteratorAtomDependent,
        AtomsetIteratorDirectable {

    /**
     * Default instance uses the atom's sequencer as the linker to begin list
     * iteration.
     */
    public AtomIteratorSequencerList() {
        this(DEFAULT);
    }

    /**
     * Constructs new class with hasNext as false. Must invoke setAtom and reset
     * before beginning iteration. Argument is class that identifies the linker
     * for beginning iteration, given an atom.
     */
    //see etomica.nbr.cell.ApiIntraspecies1ACell for non-default use of this
    // constructor
    public AtomIteratorSequencerList(AtomToLinker atomToLinker) {
        this.atomToLinker = atomToLinker;
        listIterator = new AtomIteratorList();
        setDirection(null);
    }

    public boolean hasNext() {
        return listIterator.hasNext();
    }

    public AtomSet next() {
        return nextAtom();
    }

    public Atom nextAtom() {
        Atom nextAtom = listIterator.nextAtom();
        if (doGoDown && !listIterator.hasNext())
            resetDown();
        return nextAtom;
    }

    public AtomSet peek() {
        return listIterator.peek();
    }

    public void unset() {
        listIterator.unset();
    }

    /**
     * Overrides superclass reset to ensure no reset is performed if a firstAtom
     * has not been identified. Otherwise readies for iteration beginning with
     * the numToSkip atom after firstAtom.
     */
    public void reset() {
        listIterator.unset();
        if (firstAtom == null) {
            return;
        }
        upListNow = doBoth || (direction == IteratorDirective.UP);
        doGoDown = doBoth || (direction == IteratorDirective.DOWN);

        if (upListNow) {
            resetUp();
        }
        if (!listIterator.hasNext() && doGoDown) {
            resetDown();
        }
    }

    public void allAtoms(AtomsetAction action) {
        reset();
        //can't use listIterator.allAtoms because of skipping
        while (listIterator.hasNext()) {
            action.actionPerformed(listIterator.next());
        }
        if (doGoDown) {
            resetDown();
        }
        while (listIterator.hasNext()) {
            action.actionPerformed(listIterator.next());
        }
    }

    public boolean contains(AtomSet atom) {
        if (firstAtom == null) {
            return false;
        }
        int index = listIterator.getList().indexOf((Atom) atom);
        if (index == -1) {
            return false; //not in list
        }
        int index0 = listIterator.getList().indexOf(firstAtom);
        int diff = index - index0;
        if (direction == IteratorDirective.UP) {
            return (diff - numToSkip >= 0);
        }
        if (direction == IteratorDirective.DOWN) {
            return (diff + numToSkip >= 0);
        }
        return Math.abs(diff) - numToSkip >= 0;
    }

    public int size() {
        if (firstAtom == null) {
            return 0;
        }
        int index0 = listIterator.getList().indexOf(firstAtom);
        int listSize = listIterator.getList().size();
        if (direction == IteratorDirective.UP) {
            return Math.max(0, listSize - index0 - numToSkip);
        }
        if (direction == IteratorDirective.DOWN) {
            return Math.max(0, 1 + index0 - numToSkip);
        }
        if (numToSkip == 0) {
            return listSize;
        }
        return Math.max(0, listSize - index0 - numToSkip)
                + Math.max(0, index0 - numToSkip);
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
        firstAtom = atom;
        if (atom != null) {
            listIterator.setFirst(atomToLinker.getLinker(atom));
        } else {
            listIterator.unset();
        }
    }

    /**
     * @return Returns the number of atoms to skip.
     */
    public int getNumToSkip() {
        return numToSkip;
    }

    /**
     * Sets iterator to skip some number of atoms from and including the atom
     * given by setAtom. If given 0, first atom will be the one specified in
     * setAtom. If given 1, the set-atom is skipped, and the first iterate will
     * be the atom adjacent to it up and/or down the list. Skipping is performed
     * equally in both directions, so for example if direction == null and
     * skipping is 2, the set atom and each adjacent atom will be skipped. If
     * skipping is 0 and direction == null, set-atom will be given only once (as
     * the first iterate).
     * 
     * @param numToSkip:
     *            the number of atoms to skip. Must be non-negative.
     */
    public void setNumToSkip(int numToSkip) {
        if (numToSkip < 0) {
            throw new IllegalArgumentException(
                    "Cannot give negative number to skip");
        }
        this.numToSkip = numToSkip;
    }

    /**
     * Indicates direction of iteration from target atom. A null argument
     * indicates iteration in both directions, first up from (and including, if
     * numToSkip = 0) the set-atom, then down from it (and not including it, if
     * doing both directions).
     */
    public void setDirection(IteratorDirective.Direction direction) {
        this.direction = direction;
        doBoth = (direction == null);
    }

    public IteratorDirective.Direction getDirection() {
        return direction;
    }

    /**
     * Used to set iterator for up-half of iteration only.
     */
    private void resetUp() {
        listIterator.setDirection(IteratorDirective.UP);
        listIterator.reset();
        for (int n = numToSkip; n != 0 && listIterator.hasNext(); n--) {
            listIterator.next();
        }
    }

    /**
     * Used to set iterator for down-half of iteration only.
     */
    protected void resetDown() {
        listIterator.setDirection(IteratorDirective.DOWN);
        listIterator.reset();
        for (int n = numToSkip; n != 0 && listIterator.hasNext(); n--) {
            listIterator.next();
        }

        //skip first atom if it was already given by uplist iteration
        if (upListNow && numToSkip == 0) {
            listIterator.next();
        }
        upListNow = false;
        doGoDown = false;
    }

    protected final AtomIteratorList listIterator;
    private int numToSkip = 0;
    protected boolean doBoth;
    protected boolean upListNow, doGoDown;
    private IteratorDirective.Direction direction;
    private Atom firstAtom = null;
    private final AtomToLinker atomToLinker;

    /**
     * Interface for class the determines the atom linker given an atom to begin
     * iteration.
     */
    public interface AtomToLinker {

        public AtomLinker getLinker(Atom atom);
    }

    private static final AtomToLinker DEFAULT = new AtomToLinker() {

        public AtomLinker getLinker(Atom atom) {
            return atom.seq;
        }
    };
}