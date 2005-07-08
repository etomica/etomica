package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;
import etomica.IteratorDirective;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.AtomLinker;

/**
 * Iterator for looping through the sequence list relative to a specified atom.
 * Iteration terminates at the first encounter of a null-atom linker.
 * Can be configured to loop upList, downList, or both direction from specified
 * atom. Also can be configured to skip any number of atoms from the specified
 * atom. No list specification is provided; rather iteration is performed
 * beginning with a linker associated with the atom via an AtomToLinker instance
 * given at construction. Default version uses the atom's sequencer as the
 * linker.
 */

/*
 * History Created on Aug 26, 2004 by kofke
 */
public class AtomIteratorSequenceDirectable implements
        AtomIteratorAtomDependent, AtomsetIteratorDirectable, java.io.Serializable {

    /**
     * Default instance uses the atom's sequencer as the linker to begin list
     * iteration.
     */
    public AtomIteratorSequenceDirectable() {
        upListIterator = new AtomIteratorSequence(IteratorDirective.UP);
        dnListIterator = new AtomIteratorSequence(IteratorDirective.DOWN);
        setDirection(null);
        iterator = upListIterator;
    }

    /**
     * Constructs new class with hasNext as false. Must invoke setAtom or
     * setFirst and reset before beginning iteration. AtomToLinker instance
     * provides a rule that identifies the linker for beginning iteration, given
     * an atom.
     * 
     * @throws IllegalArgumentException
     *             if direction is null
     */
    //see etomica.nbr.cell.ApiIntraspecies1ACell for non-default use of this
    // constructor
    public AtomIteratorSequenceDirectable(
            AtomIteratorSequence.AtomToLinker atomToLinker) {
        upListIterator = new AtomIteratorSequence(IteratorDirective.UP,
                atomToLinker);
        dnListIterator = new AtomIteratorSequence(IteratorDirective.DOWN,
                atomToLinker);
        setDirection(null);
        iterator = upListIterator;
    }

    public boolean hasNext() {
        return iterator.hasNext();
    }

    public final AtomSet next() {
        return nextAtom();
    }

    public Atom nextAtom() {
        Atom nextAtom = iterator.nextAtom();
        if (doGoDown && !iterator.hasNext()) {
            resetDown();
        }
        return nextAtom;
    }

    public AtomSet peek() {
        return iterator.peek();
    }

    public void unset() {
        doGoDown = false;
        iterator.unset();
    }

    /**
     * Readies for iteration beginning with the numToSkip atom after firstAtom.
     */
    public void reset() {
        upListIterator.unset();
        dnListIterator.unset();

        doGoDown = (direction != IteratorDirective.UP);//handles direction ==
                                                       // null

        if (direction != IteratorDirective.DOWN) {//upListNow, handles
                                                  // direction == null
            resetUp();
        }
        if (!upListIterator.hasNext() && doGoDown) {
            resetDown();
        }
    }

    public void allAtoms(AtomsetAction action) {
        //can't use listIterator.allAtoms because of skipping
        if (direction != IteratorDirective.DOWN) {
            resetUp();
            while (upListIterator.hasNext()) {
                action.actionPerformed(upListIterator.nextAtom());
            }
        }
        if (direction != IteratorDirective.UP) {
            resetDown();
            while (dnListIterator.hasNext()) {
                action.actionPerformed(dnListIterator.nextAtom());
            }
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
        upListIterator.setAtom(atom);
        dnListIterator.setFirst(upListIterator.getFirst());
        unset();
    }

    public void setFirst(AtomLinker first) {
        upListIterator.setFirst(first);
        dnListIterator.setFirst(first);
        unset();
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
    }

    public IteratorDirective.Direction getDirection() {
        return direction;
    }

    /**
     * Used to set iterator for up-half of iteration only.
     */
    private void resetUp() {
        upListIterator.reset();
        for (int n = numToSkip; n != 0; n--) {//ok if iterator expires before n = 0
            upListIterator.nextAtom();
        }
        iterator = upListIterator;
    }

    /**
     * Used to set iterator for down-half of iteration only.
     */
    protected void resetDown() {
        dnListIterator.reset();
        for (int n = numToSkip; n != 0; n--) {//ok if iterator expires before n = 0
            dnListIterator.nextAtom();
        }

        //skip first atom if it was already given by uplist iteration
        if ((direction != IteratorDirective.DOWN) && (numToSkip == 0)) {
            dnListIterator.nextAtom();
        }
        doGoDown = false;
        iterator = dnListIterator;
    }

    private final AtomIteratorSequence upListIterator, dnListIterator;
    protected AtomIteratorSequence iterator;
    private int numToSkip = 0;
    protected boolean doGoDown;
    private IteratorDirective.Direction direction;
    private final AtomsetCount counter = new AtomsetCount();
    private final AtomsetDetect detector = new AtomsetDetect(null);

}
