package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Pair iterator synthesized from two atom iterators, such that the inner-loop
 * iteration depends on the outer-loop atom. Pairs are formed from the atoms
 * yielded by the two atom iterators. The inner-loop iterator must implement
 * AtomIteratorAtomDependent, and its set(Atom) method will be invoked with the
 * current outer-loop atom at the start of each inner-loop iteration. All pairs
 * returned by iterator are the same Atom[] instance, and differ only in the
 * Atom instances held by it.
 * <p>
 * Iterator can be condition to put the atoms in either order in the AtomPair that
 * it returns.  Thus the inner-loop Atom may be atom0 of the returned AtomPair,
 * or it may be atom1.  This behavior is set at construction, and cannot be changed
 * afterwards.  Default behavior has outer loop atoms as atom0, and inner loop atoms
 * as atom1.
 */
public class ApiInnerVariable implements AtomPairIterator, ApiComposite, java.io.Serializable {

    /**
     * Construct a pair iterator using the given atom iterators. Requires call
     * to reset() before beginning iteration.
     */
    public ApiInnerVariable(AtomIterator aiOuter,
            AtomIteratorAtomDependent aiInner) {
        this(aiOuter,aiInner,false);
    }
    
    /**
     * Construct a pair iterator using the given atom iterators, indicating
     * whether the atoms ordering in the AtomPair should be swapped from the
     * default behavior.
     * 
     * @param aiOuter
     *            outer-loop iterator
     * @param aiInner
     *            inner-loop iterator
     * @param doSwap
     *            if false (default), outer-loop atoms are given in atom0, and
     *            inner loop in atom1; if true outer-loop atoms are given in
     *            atom1, and inner loop in atom0
     */
    public ApiInnerVariable(AtomIterator aiOuter,
            AtomIteratorAtomDependent aiInner, boolean doSwap) {
        this.aiOuter = aiOuter;
        this.aiInner = aiInner;
        this.doSwap = doSwap;
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
            if (doSwap) {
                pair.atom1 = aiOuter.nextAtom();
            }
            else {
                pair.atom0 = aiOuter.nextAtom();
            }
            aiInner.setAtom(doSwap ? pair.atom1 : pair.atom0);
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
            if (doSwap) {
                pair.atom1 = atom1;
            }
            else {
                pair.atom0 = atom1;
            }
            needUpdate1 = false;
        } //aiOuter was advanced
        if (doSwap) {
            pair.atom0 = (Atom) aiInner.peek();
        }
        else {
            pair.atom1 = (Atom) aiInner.peek();
        }
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
            if (doSwap) {
                pair.atom1 = atom1;
            }
            else {
                pair.atom0 = atom1;
            }
            needUpdate1 = false;
        } //aiOuter was advanced
        if (doSwap) {
            pair.atom0 = aiInner.nextAtom();
        }
        else {
            pair.atom1 = aiInner.nextAtom();
        }
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
            if (doSwap) {
                pair.atom1 = aiOuter.nextAtom();
                aiInner.setAtom(pair.atom1);
                aiInner.reset();
                while (aiInner.hasNext()) {
                    pair.atom0 = aiInner.nextAtom();
                    act.actionPerformed(pair);
                }
            }
            else {
                pair.atom0 = aiOuter.nextAtom();
                aiInner.setAtom(pair.atom0);
                aiInner.reset();
                while (aiInner.hasNext()) {
                    pair.atom1 = aiInner.nextAtom();
                    act.actionPerformed(pair);
                }
            }
        }
    }

    public final int nBody() {
        return 2;
    }

    private static final long serialVersionUID = 1L;
    protected final AtomPair pair = new AtomPair();
    protected boolean hasNext, needUpdate1;
    protected Atom atom1;

    /**
     * The iterators used to generate the sets of atoms.
     */
    protected final AtomIteratorAtomDependent aiInner;
    protected final AtomIterator aiOuter;
    protected final boolean doSwap;

}

