/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Phase;
import etomica.IteratorDirective.Direction;
import etomica.action.AtomsetAction;
import etomica.atom.AtomPair;
import etomica.atom.AtomPairIterator;
import etomica.atom.AtomSet;

/**
 * Adapater class that wraps three atomPair iterators, one suitable for
 * iterating over all molecule pairs in a phase (AA), another suitable for
 * iterating over all molecule pairs formed with a target molecule (1A), and the
 * third suitable for iterating over a single molecule pair (11). Appropriate
 * iterator is selected based on argument given to setTarget method. If
 * setTarget method is never called, or has count == 0, the all-pair iterator is
 * used by default. <p>
 * This class may be set up to do inter- or intra-species iteration, depending
 * on choice of inner iterators given at construction.  Wrapped iterators are final
 * and cannot be changed after construction.
 * 
 */
public class ApiMolecule implements AtomsetIteratorMolecule, AtomPairIterator, java.io.Serializable {

    /**
     * Constructs iterator by wrapping three others.
     * @param api11
     *            iterator for single molecule pair formed from target
     *            specification
     * @param api1A
     *            iterator for all pairs formed with a target molecule
     * @param apiAA
     *            iterator for all pairs in the phase
     */
    public ApiMolecule(AtomsetIteratorMolecule api11,
            AtomsetIteratorMolecule api1A, AtomsetIteratorPhaseDependent apiAA) {
        this.api11 = api11;
        this.api1A = api1A;
        this.apiAA = apiAA;
        iterator = (AtomPairIterator) apiAA;
    }

    /**
     * Sets target atoms and determines the pair iterator according to the value
     * of targetAtoms.count(). Thus, for count equal to
     * <ul>
     * <li>0, AA iterator is indicated
     * <li>1, 1A iterator is indicated
     * <li>2 or larger, 11 iterator is indicated
     * </ul>
     * Target is passed on to setTarget method of corresponding iterator, which
     * ultimately determines the behavior of this iterator.
     * 
     * @throws NullPointerException
     *             if targetAtoms is null; use AtomSet.NULL instead
     */
    public void setTarget(AtomSet targetAtoms) {
        switch (targetAtoms.count()) {
        case 0:
            iterator = (AtomPairIterator) apiAA;
            break;
        case 1:
            iterator = (AtomPairIterator) api1A;
            api1A.setTarget(targetAtoms);
            break;
        default:
            iterator = (AtomPairIterator) api11;
            api11.setTarget(targetAtoms);
        }
        getCurrentIterator().setPhase(phase);
    }

    /**
     * Specifies the phase from which iterates are taken.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        getCurrentIterator().setPhase(phase);
    }

    /**
     * Specifies the direction from which partners of a target atoms
     * are taken to form iterates.  Has no effect if performing AA or 11
     * iteration (that is, if set target has count not equal to 1).
     */
    public void setDirection(Direction direction) {
        api1A.setDirection(direction);
    }

    /**
     * Returns true if the iterator as currently conditioned will return
     * the given set of atoms.  Does not require reset; clobbers iteration
     * state.
     */
    public boolean contains(AtomSet atom) {
        return iterator.contains(atom);
    }

    /**
     * Indicates if iterator has another iterate.
     */
    public boolean hasNext() {
        return iterator.hasNext();
    }

    /**
     * Readies iterator for iteration.
     */
    public void reset() {
        iterator.reset();
    }

    /**
     * Sets iterator to state in which hasNext is false.
     */
    public void unset() {
        iterator.unset();
    }

    /**
     * Same as nextPair.
     */
    public AtomSet next() {
        return nextPair();
    }

    /**
     * Returns next iterate.
     */
    public AtomPair nextPair() {
        return iterator.nextPair();
    }

    /**
     * Returns next iterate without advancing iterator.
     */
    public AtomSet peek() {
        return iterator.peek();
    }

    /**
     * Performs action on all iterates.  Does not require
     * reset; clobbers iteration state.
     */
    public void allAtoms(AtomsetAction action) {
        iterator.allAtoms(action);
    }

    /**
     * Returns the number of iterates that would be given
     * when iterating after reset.  Does not require reset;
     * clobbers iteration state.
     */
    public int size() {
        return iterator.size();
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public final int nBody() {
        return 2;
    }

    /**
     * Returns the iterator currently being used for the pair iteration. If a
     * target had been specified, this will be the Api1A iterator, otherwise it
     * will be the ApiAA iterator.
     */
    public AtomsetIteratorPhaseDependent getCurrentIterator() {
        return (AtomsetIteratorPhaseDependent) iterator;
    }

    /**
     * Returns the 11 iterator set at construction.
     */
    public AtomsetIteratorMolecule getApi11() {
        return api11;
    }

    /**
     * Returns the 1A iterator set at construction.
     */
    public AtomsetIteratorMolecule getApi1A() {
        return api1A;
    }

    /**
     * Returns the AA iterator set at construction.
     */
    public AtomsetIteratorPhaseDependent getApiAA() {
        return apiAA;
    }
    
    private AtomPairIterator iterator;
    private final AtomsetIteratorMolecule api11;
    private final AtomsetIteratorMolecule api1A;
    private final AtomsetIteratorPhaseDependent apiAA;
    private Phase phase;

}
