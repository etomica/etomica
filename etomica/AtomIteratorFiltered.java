/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
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
		this(iterator, AtomFilter.ALL);
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
	public boolean contains(Atom atom) {
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
	public Atom next() {
		Atom nextAtom = next;
		next = null;
		while(iterator.hasNext() && next == null) {
			next = iterator.next();
			if(!filter.accept(next)) next = null;
		}
		return nextAtom;
	}
	
	/**
	 * Returns next atom without advancing the iterator.
	 */
	public Atom peek() {
		return next;
	}

	/**
	 * Performs the given action on all atoms from iterator that 
	 * meet the filter criteria.
	 */
	public void allAtoms(AtomActive action) {
		iterator.allAtoms(actionWrapper(filter, action));
	}

	/**
	 * Returns the number of iterates given by the
	 * iterator that meet the criteria of the 
	 * filter.
	 */
	public int size() {
		AtomActiveCount counter = new AtomActiveCount();
		allAtoms(counter);
		return counter.callCount();
	}
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

	/**
	 * Returns a new action that wraps the given action such that action is performed
	 * only on the atoms meeting the filter's criteria.
	 */
	private static AtomActive actionWrapper(final AtomFilter filter, final AtomActive action) {
		return new AtomActive() {
			public void actionPerformed(Atom atom) {
				if(filter.accept(atom)) action.actionPerformed(atom);
			}
		};
	}
}
