/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomSet;
import etomica.AtomsetIterator;
import etomica.IteratorDirective.Direction;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.action.AtomsetCount;
import etomica.atom.AtomPairVector;
import etomica.atom.AtomsetFilter;

/**
 * Wraps an AtomIterator and filters its iterates so that
 * only those meeting specified criteria are returned.
 */

public class AtomsetIteratorFiltered implements AtomsetIteratorTargetable, AtomsetIteratorDirectable, 
            AtomsetIteratorBasisDependent, AtomPairIterator {
	/**
	 * Default constructor that causes no atoms to be filtered.
	 * Iterator will give all iterates of the given iterator
	 * until another filter is specified.
	 */
	public AtomsetIteratorFiltered(AtomsetIterator iterator) {
		this(iterator, AtomsetFilter.ACCEPT_ALL);
	}
	
	/**
	 * Returns the iterates of the given iterator that meet
	 * the critertia of the given filter.
	 * @param iterator
	 * @param filter
	 */
	public AtomsetIteratorFiltered(AtomsetIterator iterator, AtomsetFilter filter) {
		this.iterator = iterator;
		this.filter = filter;
        nBody = iterator.nBody();
        if (nBody == 2) {
            nextAtoms = new AtomPairVector();
            next = new AtomPairVector();
        }
        else {
            nextAtoms = null;
        }
	}

	
	/**
	 * Returns true if the iterator contains the given atom and
	 * atom meets the filter's criteria.
	 */
	public boolean contains(AtomSet atom) {
		return filter.accept(atom) && iterator.contains(atom);
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
            if (nBody == 1) {
                next = iterator.next();
            }
            else {
                ((AtomPairIterator)iterator).nextPair().copyTo((AtomPair)next);
            }
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
        if (nBody == 1) {
            nextAtom = (Atom)next;
        }
        else {
            ((AtomPair)next).copyTo(nextAtoms);
        }
        hasNext = false;
		while(iterator.hasNext()) {
            if (nBody == 1) {
                next = iterator.next();
            }
            else {
                ((AtomPairIterator)iterator).nextPair().copyTo((AtomPair)next);
            }
            if(filter.accept(next)) {
                hasNext = true;
                break;
            }
		}
        if (nBody == 1) {
            return nextAtom;
        }
        else {
            return nextAtoms;
        }
	}
	
    public AtomPair nextPair() {
        ((AtomPair)next).copyTo(nextAtoms);
        hasNext = false;
        while(iterator.hasNext()) {
            ((AtomPairIterator)iterator).nextPair().copyTo((AtomPair)next);
            if(filter.accept(next)) {
                hasNext = true;
                break;
            }
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
	
	public int basisSize() {
        return ((AtomsetIteratorBasisDependent)iterator).basisSize();
    }
    
    public void setBasis(AtomSet atoms) {
        ((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);
        if (atoms instanceof AtomPairVector) {
            ((AtomPairVector)next).nearestImageVector = ((AtomPairVector)atoms).nearestImageVector;
        }
    }
        
	/**
	 * @return Returns the filter.
	 */
	public AtomsetFilter getFilter() {
		return filter;
	}
	/**
	 * @param filter The filter to set.
	 */
	public void setFilter(AtomsetFilter filter) {
		this.filter = filter;
	}
	
	/**
	 * @return the iterator wrapped by this filter.
	 */
	public AtomsetIterator getIterator() {
		return iterator;
	}
	
	public void setDirection(Direction direction) {
		if (iterator instanceof AtomsetIteratorDirectable) {
		    ((AtomsetIteratorDirectable)iterator).setDirection(direction);
        }
	}
	public void setTarget(AtomSet targetAtoms) {
        if (iterator instanceof AtomsetIteratorTargetable) {
            ((AtomsetIteratorTargetable)iterator).setTarget(targetAtoms);
        }
	}

	private final AtomsetIterator iterator;
	private AtomsetFilter filter;
	private AtomSet next;
	private final AtomPairVector nextAtoms;
    private final int nBody;
    private Atom nextAtom;
    private boolean hasNext;

	/**
	 * Returns a new action that wraps the given action such that action is performed
	 * only on the atoms meeting the filter's criteria.
	 */
	private static AtomsetAction actionWrapper(final AtomsetFilter filter, final AtomsetAction action) {
		return new AtomsetActionAdapter() {
			public void actionPerformed(AtomSet atom) {
				if(filter.accept(atom)) action.actionPerformed(atom);
			}
		};
	}
}
