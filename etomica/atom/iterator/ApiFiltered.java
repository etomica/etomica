package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.action.AtomsetCount;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.atom.AtomsetFilter;
import etomica.atom.iterator.IteratorDirective.Direction;

/**
 * Wraps an AtomPairIterator and filters its iterates so that
 * only those meeting specified criteria are returned.
 */

public class ApiFiltered implements AtomsetIteratorDirectable, 
            AtomsetIteratorBasisDependent, java.io.Serializable {
	
    /**
	 * Returns the iterates of the given iterator that meet
	 * the critertia of the given filter.
	 */
	public ApiFiltered(AtomsetIterator iterator, AtomsetFilter filter) {
		this.iterator = iterator;
		this.filter = filter;
	}

    /**
     * Puts iterator in state ready for iteration.
     */
    public void reset() {
	    iterator.reset();
    }

	/**
	 * Sets iterator so that hasNext returns false.
	 */
	public void unset() {
		iterator.unset();
	}

	/**
	 * Returns the next atom from the iterator that meets the 
	 * filter's criteria.
	 */
	public IAtomSet next() {
        IAtomSet next = iterator.next();
        while (next != null && !filter.accept(next)) {
            next = iterator.next();
        }
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
    public void setBasis(IAtomSet atoms) {
        ((AtomsetIteratorBasisDependent)iterator).setBasis(atoms);
    }
    
    /**
     * Passes call on to wrapped iterator and returns value.
     */
    public boolean haveTarget(IAtom target) {
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
	public void setTarget(IAtom targetAtom) {
        if (iterator instanceof AtomsetIteratorTargetable) {
            ((AtomsetIteratorTargetable)iterator).setTarget(targetAtom);
        }
	}

    private static final long serialVersionUID = 2L;
	protected final AtomsetIterator iterator;
	protected final AtomsetFilter filter;

	/**
	 * Returns a new action that wraps the given action such that action is performed
	 * only on the atoms meeting the filter's criteria.
	 */
	private static AtomsetAction actionWrapper(final AtomsetFilter filter, final AtomsetAction action) {
		return new AtomsetActionAdapter() {
			public void actionPerformed(IAtomSet atom) {
				if(filter.accept(atom)) action.actionPerformed(atom);
			}
		};
	}
}
