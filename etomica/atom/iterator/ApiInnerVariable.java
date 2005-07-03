package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomSet;
import etomica.action.AtomsetAction;

/**
 * Pair iterator synthesized from two atom iterators, such that the inner-loop
 * iteration depends on the outer-loop atom. Pairs are formed from the atoms
 * yielded by the two atom iterators. The inner-loop iterator must implement
 * AtomIteratorAtomDependent, and its set(Atom) method will be invoked with the
 * current outer-loop atom at the start of each inner-loop iteration. All pairs
 * returned by iterator are the same Atom[] instance, and differ only in the
 * Atom instances held by it.
 */

/*
 * History of changes 08/25/04 (DAK et al) new
 */

public class ApiInnerVariable implements AtomPairIterator, ApiComposite {

    /**
     * Construct a pair iterator using the given atom iterators. Requires call
     * to reset() before beginning iteration.
     */
    public ApiInnerVariable(AtomIterator aiOuter,
            AtomIteratorAtomDependent aiInner) {
        this.aiOuter = aiOuter;
        this.aiInner = aiInner;
        unset();
    }

    /**
     * Accessor method for the outer-loop atom iterator.
     * 
     * @return the current outer-loop iterator
     */
    public AtomIterator getOuterIterator() {
        return aiOuter;
    }

    /**
     * Accessor method for the inner-loop atom iterator.
     * 
     * @return the current inner-loop iterator
     */
    public AtomIterator getInnerIterator() {
        return aiInner;
    }

    /**
     * Sets the iterator such that hasNext is false.
     */
    public void unset() {
        hasNext = false;
    }

    /**
     * Indicates whether the given atom pair will be returned by the iterator
     * during its iteration. The order of the atoms in the pair is significant
     * (this means that a value of true is returned only if one of the pairs
     * returned by the iterator will have the same two atoms in the same
     * atom1/atom2 position as the given pair). Not dependent on state of
     * hasNext.
     */
    public boolean contains(AtomSet pair) {
        if(pair == null || pair.count() != 2) {
            return false;
        }
        if (aiOuter.contains(pair.getAtom(0))) {
            aiInner.setAtom(pair.getAtom(0));
            return aiInner.contains(pair.getAtom(1));
        }
        return false;
    }

    /**
     * Returns the number of pairs given by this iterator. Independent of state
     * of hasNext. Clobbers the iteration state (i.e., status of hasNext/next)
     * but does not recondition iterator (i.e., does not change set of iterates
     * that would be given on iteration after reset). Must perform reset if
     * attempting iteration after using size() method.
     */
    public int size() {
        int sum = 0;
        aiOuter.reset();
        while (aiOuter.hasNext()) {
            aiInner.setAtom(aiOuter.nextAtom());
            sum += aiInner.size();
        }
        return sum;
    }

    /**
     * Indicates whether the iterator has completed its iteration.
     */
    public boolean hasNext() {
        return hasNext;
    }

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        aiOuter.reset();
        hasNext = false;
        needUpdate1 = false;
        while (aiOuter.hasNext()) { //loop over outer iterator...
            pair.atom0 = aiOuter.nextAtom();
            aiInner.setAtom(pair.atom0);
            aiInner.reset();
            if (aiInner.hasNext()) { //until inner iterator has another
                hasNext = true;
                break; //...until iterator 2 hasNext
            }
        }//end while
    }

    /**
     * Returns the next pair without advancing the iterator. If the iterator has
     * reached the end of its iteration, returns null.
     */
    public AtomSet peek() {
        if (!hasNext) {
            return null;
        }
        if (needUpdate1) {
            pair.atom0 = atom1;
            needUpdate1 = false;
        } //aiOuter was advanced
        pair.atom1 = (Atom) aiInner.peek();
        return pair;
    }

    public AtomSet next() {
        return nextPair();
    }

    /**
     * Returns the next pair of atoms. The same Atom[] instance is returned
     * every time, but the Atoms it holds are (of course) different for each
     * iterate.
     */
    public AtomPair nextPair() {
        if (!hasNext)
            return null;
        //we use this update flag to indicate that atom1 in pair needs to be
        // set to a new value.
        //it is not done directly in the while-loop because pair must first
        // return with the old atom1 intact
        if (needUpdate1) {
            pair.atom0 = atom1;
            needUpdate1 = false;
        } //aiOuter was advanced
        pair.atom1 = aiInner.nextAtom();
        while (!aiInner.hasNext()) { //Inner is done for this atom1, loop until
                                     // it is prepared for next
            if (aiOuter.hasNext()) { //Outer has another atom1...
                atom1 = aiOuter.nextAtom(); //...get it
                aiInner.setAtom(atom1);
                aiInner.reset();
                needUpdate1 = true; //...flag update of pair.atom1 for next
                                    // time
            } else {
                hasNext = false;
                break;
            } //Outer has no more; all done with pairs
        }//end while
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allAtoms(AtomsetAction act) {
        aiOuter.reset();
        while (aiOuter.hasNext()) {
            pair.atom0 = aiOuter.nextAtom();
            aiInner.setAtom(pair.atom0);
            aiInner.reset();
            while (aiInner.hasNext()) {
                pair.atom1 = aiInner.nextAtom();
                act.actionPerformed(pair);
            }
        }
    }

    public final int nBody() {
        return 2;
    }

    protected final AtomPair pair = new AtomPair();
    protected boolean hasNext, needUpdate1;
    protected Atom atom1;

    /**
     * The iterators used to generate the sets of atoms.
     */
    protected final AtomIteratorAtomDependent aiInner;
    protected final AtomIterator aiOuter;

} //end of class ApiInnerVariable

