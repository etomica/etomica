/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomSet;
import etomica.AtomsetIterator;
import etomica.IteratorDirective.Direction;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.action.AtomsetCount;
import etomica.atom.AtomsetFilter;

/**
 * Wraps an AtomIterator and filters its iterates so that
 * only those meeting specified criteria are returned.
 */

//XXX this shouldn't have to be AtomsetIteratorBasisDependent, but as used it does
public class AtomsetIteratorFiltered implements AtomsetIteratorTargetable, AtomsetIteratorDirectable, 
            AtomsetIteratorBasisDependent {
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
            nextAtoms = new AtomPair();
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
		return next != null;
	}

	/**
	 * Puts iterator in state ready for iteration.
	 */
	public void reset() {
		iterator.reset();
		next = null;
		while(iterator.hasNext() && next == null) {
			next = iterator.next();
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

	/**
	 * Returns the next atom from the iterator that meets the 
	 * filter's criteria.
	 */
	public AtomSet next() {
        if (nBody == 1) {
            nextAtom = (Atom)next;
        }
        else {
            ((AtomPair)next).copyTo((AtomPair)nextAtoms);
        }
		next = null;
		while(iterator.hasNext() && next == null) {
			next = iterator.next();
			if(!filter.accept(next)) next = null;
		}
        if (nBody == 1) {
            return nextAtom;
        }
        else {
            return nextAtoms;
        }
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
	private final AtomSet nextAtoms;
    private final int nBody;
    private Atom nextAtom;

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
