/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica.nbr;

import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomSet;
import etomica.Phase;
import etomica.IteratorDirective.Direction;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.action.AtomsetCount;
import etomica.atom.AtomPairFilter;
import etomica.atom.AtomPairVector;
import etomica.atom.iterator.ApiMolecule;
import etomica.atom.iterator.AtomsetIteratorMolecule;

/**
 * Wraps an AtomsetIterator and filters its iterates so that
 * only those meeting specified criteria are returned.  Explicitly
 * written for case where atom set is a pair.
 */
//maybe make an AtomsetIteratorMoleculeFiltered that implements AtomsetIteratorMolecule
//and leave this to implement just AtomsetIterator
public class ApiFiltered implements AtomsetIteratorMolecule, AtomPairIterator {
    /**
     * Default constructor that causes no atoms to be filtered.
     * Iterator will give all iterates of the given iterator
     * until another filter is specified.
     */
    public ApiFiltered(ApiMolecule iterator) {
        this(iterator, new NeighborCriterionAll());
    }

    /**
     * Returns the iterates of the given iterator that meet
     * the critertia of the given filter.
     * @param iterator
     * @param filter
     */
    public ApiFiltered(ApiMolecule iterator, NeighborCriterion filter) {
        if(iterator.nBody() != 2) throw new IllegalArgumentException("Illegal attempt to construct pair iterator by wrapping a non-pair iterator");
        this.iterator = iterator;
        this.filter = filter;
        nextAtoms = new AtomPairVector();
    }


    /**
     * Returns true if the iterator contains the given atom and
     * atom meets the filter's criteria.
     */
    public boolean contains(AtomSet atom) {
        return filter.accept((AtomPair)atom) && iterator.contains(atom);
    }

    /**
     * Indicates whether iterator has another iterate to return.
     */
    public boolean hasNext() {
        return next != null;
    }

    /**
     * Puts iterator in state ready for iteration.
     */
    public void reset() {
        iterator.reset();
        next = null;
        while(iterator.hasNext() && next == null) {
            next = (AtomPairVector)iterator.next();
            if(!filter.accept(next)) next = null;
        }
    }
    
    /**
     * Sets iterator so that hasNext returns false.
     */
    public void unset() {
        iterator.unset();
        next = null;
    }
    
    public AtomSet next() {
        return nextPair();
    }

    /**
     * Returns the next atom from the iterator that meets the 
     * filter's criteria.
     */
    public AtomPair nextPair() {
        if(next == null) return null;
        next.copyTo(nextAtoms);
        next = null;
        while(iterator.hasNext() && next == null) {
            next = iterator.nextPair();
            if(!filter.accept(next)) next = null;
        }
        return nextAtoms;
    }

    /**
     * Returns next atom without advancing the iterator.
     */
    public AtomSet peek() {
        return next;
    }

    /**
     * Performs the given action on all atoms from iterator that 
     * meet the filter criteria.
     */
    public void allAtoms(AtomsetAction action) {
        iterator.allAtoms(actionWrapper(filter, action));
    }

    /**
     * Returns the number of iterates given by the
     * iterator that meet the criteria of the 
     * filter.
     */
    public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
    }

    public final int nBody() {
        return iterator.nBody();
    }

    /**
     * @return Returns the filter.
     */
    public AtomPairFilter getFilter() {
        return filter;
    }

    /**
     * @return the iterator wrapped by this filter.
     */
    public ApiMolecule getIterator() {
        return iterator;
    }

    public void setPhase(Phase phase) {
        iterator.setPhase(phase);
    }
    public void setDirection(Direction direction) {
        iterator.setDirection(direction);
    }
    public void setTarget(AtomSet targetAtoms) {
        iterator.setTarget(targetAtoms);
    }

    private final ApiMolecule iterator;
    private final NeighborCriterion filter;
    private AtomPair next;
    private final AtomPairVector nextAtoms;

    /**
     * Returns a new action that wraps the given action such that action is performed
     * only on the atoms meeting the filter's criteria.
     */
    private static AtomsetAction actionWrapper(final AtomPairFilter filter, final AtomsetAction action) {
        return new AtomsetActionAdapter() {
            public void actionPerformed(AtomSet atom) {
                if(filter.accept((AtomPair)atom)) action.actionPerformed(atom);
            }
        };
    }
}
