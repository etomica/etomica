/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica.nbr;

import etomica.ApiMolecule;
import etomica.Atom;
import etomica.AtomsetFilter;
import etomica.AtomsetIterator;
import etomica.AtomsetIteratorMolecule;
import etomica.Phase;
import etomica.IteratorDirective.Direction;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.action.AtomsetCount;
import etomica.lattice.RectangularLattice;
import etomica.nbr.cell.AtomsetIteratorCellular;

/**
 * Wraps an AtomsetIterator and filters its iterates so that
 * only those meeting specified criteria are returned.  Explicitly
 * written for case where atom set is a pair.
 */
//maybe make an AtomsetIteratorMoleculeFiltered that implements AtomsetIteratorMolecule
//and leave this to implement just AtomsetIterator
public class ApiFiltered implements AtomsetIteratorMolecule {
	/**
	 * Default constructor that causes no atoms to be filtered.
	 * Iterator will give all iterates of the given iterator
	 * until another filter is specified.
	 */
    //specific to AtomsetIteratorMolecule because used by PotentialMasterCell
//	public ApiFiltered(AtomsetIteratorMolecule iterator) {
//		this(iterator, AtomsetFilter.ACCEPT_ALL);
//	}
	
	/**
	 * Returns the iterates of the given iterator that meet
	 * the critertia of the given filter.
	 * @param iterator
	 * @param filter
	 */
	public ApiFiltered(ApiMolecule iterator, NeighborCriterionSimple filter) {
        if(iterator.nBody() != 2) throw new IllegalArgumentException("Illegal attempt to construct pair iterator by wrapping a non-pair iterator");
		this.iterator = iterator;
		this.filter = filter;
		nextAtoms = new Atom[2];
	}

	
	/**
	 * Returns true if the iterator contains the given atom and
	 * atom meets the filter's criteria.
	 */
	public boolean contains(Atom[] atom) {
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
        filter.setCellIterator((AtomsetIteratorCellular)iterator.getCurrentIterator());
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
	public Atom[] next() {
        if(next == null) return null;
        nextAtoms[0] = next[0];
        nextAtoms[1] = next[1];
		next = null;
		while(iterator.hasNext() && next == null) {
			next = iterator.next();
			if(!filter.accept(next)) next = null;
		}
		return nextAtoms;
	}
	
	/**
	 * Returns next atom without advancing the iterator.
	 */
	public Atom[] peek() {
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
	public AtomsetFilter getFilter() {
		return filter;
	}
	
	/**
	 * @return the iterator wrapped by this filter.
	 */
	public AtomsetIterator getIterator() {
		return iterator;
	}
	
	public void setPhase(Phase phase) {
		iterator.setPhase(phase);
	}
	public void setDirection(Direction direction) {
		iterator.setDirection(direction);
	}
	public void setTarget(Atom[] targetAtoms) {
		iterator.setTarget(targetAtoms);
	}

	private final ApiMolecule iterator;
	private final NeighborCriterionSimple filter;
	private Atom[] next;
	private final Atom[] nextAtoms;

	/**
	 * Returns a new action that wraps the given action such that action is performed
	 * only on the atoms meeting the filter's criteria.
	 */
	private static AtomsetAction actionWrapper(final AtomsetFilter filter, final AtomsetAction action) {
		return new AtomsetActionAdapter() {
			public void actionPerformed(Atom[] atom) {
				if(filter.accept(atom)) action.actionPerformed(atom);
			}
		};
	}
}
