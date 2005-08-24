/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.action.AtomsetCount;
import etomica.atom.AtomPair;
import etomica.atom.AtomPairFilter;
import etomica.atom.AtomPairVector;
import etomica.atom.AtomSet;
import etomica.atom.iterator.IteratorDirective.Direction;

/**
 * Wraps an AtomPairIterator and filters its iterates so that
 * only those meeting specified criteria are returned.
 */

public class ApiFiltered implements AtomsetIteratorDirectable, 
            AtomsetIteratorBasisDependent, AtomPairIterator, java.io.Serializable {
	
	/**
	 * Returns the iterates of the given iterator that meet
	 * the critertia of the given filter.
	 */
	public ApiFiltered(AtomPairIterator iterator, AtomPairFilter filter) {
		this.iterator = iterator;
		this.filter = filter;
        nextAtoms = new AtomPairVector();
        next = new AtomPairVector();
	}

	
	/**
	 * Returns true if the iterator contains the given atom and
	 * atom meets the filter's criteria.  Returns false if the 
     * argument is null, or is not an instance of AtomPair.
	 */
	public boolean contains(AtomSet pair) {
        if(!(pair instanceof AtomPair)) {
            return false;
        }
		return filter.accept((AtomPair)pair) && iterator.contains(pair);
	}
	
	/**
	 * Indicates whether iterator has another iterate to return.
	 */
	public boolean hasNext() {
		return hasNext;
	}

	/**
	 * Puts iterator in state ready for iteration.
	 */
	public void reset() {
		iterator.reset();
        hasNext = false;
		while(iterator.hasNext()) {
            iterator.nextPair().copyTo(next);
            if(filter.accept(next)) {
                hasNext = true;
                break;
            }
		}
	}

	/**
	 * Sets iterator so that hasNext returns false.
	 */
	public void unset() {
		iterator.unset();
		hasNext = false;
	}

	/**
	 * Returns the next atom from the iterator that meets the 
	 * filter's criteria.
	 */
	public AtomSet next() {
        next.copyTo(nextAtoms);
        hasNext = false;
		while(iterator.hasNext()) {
            iterator.nextPair().copyTo(next);
            if(filter.accept(next)) {
                hasNext = true;
                break;
            }
		}
        return nextAtoms;
	}
	
    /**
     * Returns the next pair and advances the iterator.
     */
    public AtomPair nextPair() {
        next.copyTo(nextAtoms);
        hasNext = false;
        while(iterator.hasNext()) {
            iterator.nextPair().copyTo(next);
            if(filter.accept(next)) {
                hasNext = true;
                break;
            }
        }
        return nextAtoms;
    }
    
	/**
	 * Returns next pair without advancing the iterator.
	 */
	public AtomSet peek() {
		return next;
	}

	/**
	 * Performs the given action on all pairs from iterator that 
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
	
    /**
     * Returns 2, indicating that this is a pair iterator.
     */
	public final int nBody() {
	    return 2;
	}
	
    /**
     * Passes call on to wrapped iterator and returns value.
     */
	public int basisSize() {
        return ((AtomsetIteratorBasisDependent)iterator).basisSize();
    }
    
    /**
     * Passes call on to wrapped iterator, if it implements AtomsetIteratorBasisDependent.
     */
    public void setBasis(AtomSet atoms) {
        ((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);
        if (atoms instanceof AtomPairVector) {
            ((AtomPairVector)next).nearestImageVector = ((AtomPairVector)atoms).nearestImageVector;
        }
    }
    
    /**
     * Passes call on to wrapped iterator and returns value.
     */
    public boolean haveTarget(AtomSet target) {
        return ((AtomsetIteratorBasisDependent)iterator).haveTarget(target);
    }
        	
	/**
	 * @return the iterator wrapped by this filter.
	 */
	public AtomsetIterator getIterator() {
		return iterator;
	}
	
    /**
     * Passes call on to wrapped iterator, if it implements AtomsetIteratorDirectable.
     */
	public void setDirection(Direction direction) {
		if (iterator instanceof AtomsetIteratorDirectable) {
		    ((AtomsetIteratorDirectable)iterator).setDirection(direction);
        }
	}
    
    /**
     * Passes call on to wrapped iterator, if it implements AtomsetIteratorTargetable.
     */
	public void setTarget(AtomSet targetAtoms) {
        if (iterator instanceof AtomsetIteratorTargetable) {
            ((AtomsetIteratorTargetable)iterator).setTarget(targetAtoms);
        }
	}

	private final AtomPairIterator iterator;
	private final AtomPairFilter filter;
	private AtomPair next;
	private final AtomPairVector nextAtoms;
    private boolean hasNext;

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
