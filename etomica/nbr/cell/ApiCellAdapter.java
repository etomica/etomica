/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Atom;
import etomica.AtomsetIteratorPhaseDependent;
import etomica.AtomsetIteratorTargetable;
import etomica.Phase;
import etomica.action.AtomsetAction;

/**
 * Adapater class that wraps two atomPair iterators, one suitable for
 * iterating over all molecule pairs in a phase, and the other suitable for iterating
 * over all molecule pairs formed with a target molecule.  Appropriate iterator
 * is selected based on argument given to setTarget method.  If method is
 * not called, the all-pair iterator is used by default.<br>
 * This class may be set up to do inter- or intra-species iteration, depending
 * on choice of inner iterators given at construction.
 */
public class ApiCellAdapter implements AtomsetIteratorTargetable,
                                                 AtomsetIteratorPhaseDependent {

    /**
     * @param api1A iterator for all pairs formed with a target molecule
     * @param apiAA iterator for all pairs in the phase
     */
    public ApiCellAdapter(AtomsetIteratorPhaseDependent api1A,
                          AtomsetIteratorPhaseDependent apiAA) {
		this.api1A = api1A;
        this.apiAA = apiAA;
        iterator = apiAA;
	}
    
    

    public void setTarget(Atom[] targetAtoms) {
        if(targetAtoms[0] == null) {
            iterator = apiAA;
        } else {
            iterator = api1A;
            ((AtomsetIteratorTargetable)api1A).setTarget(targetAtoms);
        }
        iterator.setPhase(phase);
    }
    
    public void setPhase(Phase phase) {
        this.phase = phase;
        iterator.setPhase(phase);
    }
    
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#contains(etomica.Atom)
	 */
	public boolean contains(Atom[] atom) {
		return iterator.contains(atom);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#hasNext()
	 */
	public boolean hasNext() {
		return iterator.hasNext();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#reset()
	 */
	public void reset() {
		iterator.reset();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#unset()
	 */
	public void unset() {
		iterator.unset();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#next()
	 */
	public Atom[] next() {
		return iterator.next();
	}
	
	public Atom[] peek() {
		return iterator.peek();
	}
	
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#all(etomica.Atom, etomica.IteratorDirective, etomica.AtomActive)
	 */
	public void allAtoms(AtomsetAction action) {
		iterator.allAtoms(action);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#size()
	 */
	public int size() {
		return iterator.size();
	}
	
	public final int nBody() {
		return 2;
	}
    
    private AtomsetIteratorPhaseDependent iterator;
    private final AtomsetIteratorPhaseDependent api1A, apiAA;
    private Phase phase;

}
