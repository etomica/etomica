/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomSet;
import etomica.action.AtomAction;
import etomica.action.AtomActionAdapter;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.atom.AtomFilter;

/**
 * Wraps an AtomIterator and filters its iterates so that
 * only those meeting specified criteria are returned.
 */
public class AtomIteratorFiltered implements AtomIterator {

	/**
	 * Returns the iterates of the given iterator that meet
	 * the critertia of the given filter.  Filter is final and
     * cannot be changed after construction.
	 * @param iterator the wrapped iterator
	 * @param filter iterator returns only those atoms for which
     * this filter's accept method returns true
	 */
	public AtomIteratorFiltered(AtomIterator iterator, AtomFilter filter) {
		this.iterator = iterator;
        this.filter = filter;
	}

	
	/**
	 * Returns true if the iterator contains the given atom and
	 * atom meets the filter's criteria.
	 */
	public boolean contains(AtomSet atom) {
		return filter.accept((Atom)atom) && iterator.contains(atom);
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
		nextAtom();
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
	public Atom nextAtom() {
		Atom nextAtom = next;
		next = null;
		while(iterator.hasNext() && next == null) {
			next = iterator.nextAtom();
			if(!filter.accept(next)) next = null;
		}
		return nextAtom;
	}
	
	public AtomSet next() {
		return nextAtom();
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
		iterator.allAtoms(actionWrapper(filter,(AtomAction)action));
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
	
	public final int nBody() {return 1;}
		
	private final AtomIterator iterator;
	private final AtomFilter filter;
	private Atom next;

	/**
	 * Returns a new action that wraps the given action such that action is performed
	 * only on the atoms meeting the filter's criteria.
	 */
	private static AtomAction actionWrapper(final AtomFilter filter, final AtomAction action) {
		return new AtomActionAdapter() {
			public void actionPerformed(Atom atom) {
				if(filter.accept(atom)) action.actionPerformed(atom);
			}
		};
	}
}
