/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.action.AtomsetCount;

/**
 * Wraps an AtomIterator and filters its iterates so that
 * only those meeting specified criteria are returned.
 */
public class AtomIteratorFiltered implements AtomIterator {

	/**
	 * Default constructor that causes no atoms to be filtered.
	 * Iterator will give all iterates of the given iterator
	 * until another filter is specified.
	 */
	public AtomIteratorFiltered(AtomIterator iterator) {
		this(iterator, AtomFilter.ACCEPT_ALL);
	}
	
	/**
	 * Returns the iterates of the given iterator that meet
	 * the critertia of the given filter.
	 * @param iterator
	 * @param filter
	 */
	public AtomIteratorFiltered(AtomIterator iterator, AtomFilter filter) {
		this.iterator = iterator;
	}

	
	/**
	 * Returns true if the iterator contains the given atom and
	 * atom meets the filter's criteria.
	 */
	public boolean contains(Atom[] atom) {
		return filter.accept(atom[0]) && iterator.contains(atom);
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
		next();
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
	
	public Atom[] next() {
		atoms[0] = nextAtom();
		return atoms;
	}
	
	/**
	 * Returns next atom without advancing the iterator.
	 */
	public Atom[] peek() {
		atoms[0] = next;
		return atoms;
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
	
	public final int nBody() {return 1;}
	
	/**
	 * @return Returns the filter.
	 */
	public AtomFilter getFilter() {
		return filter;
	}
	/**
	 * @param filter The filter to set.
	 */
	public void setFilter(AtomFilter filter) {
		this.filter = filter;
	}
	
	private final AtomIterator iterator;
	private AtomFilter filter;
	private Atom next;
	private final Atom[] atoms = new Atom[1];

	/**
	 * Returns a new action that wraps the given action such that action is performed
	 * only on the atoms meeting the filter's criteria.
	 */
	private static AtomsetAction actionWrapper(final AtomFilter filter, final AtomsetAction action) {
		return new AtomsetActionAdapter() {
			public void actionPerformed(Atom[] atom) {
				if(filter.accept(atom[0])) action.actionPerformed(atom);
			}
		};
	}
}
