/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.AtomSet;
import etomica.NearestImageVectorSource;
import etomica.Phase;
import etomica.IteratorDirective.Direction;
import etomica.action.AtomsetAction;
import etomica.space.Vector;

/**
 * Adapater class that wraps two atomPair iterators, one suitable for
 * iterating over all molecule pairs in a phase (AA), and the other suitable for iterating
 * over all molecule pairs formed with a target molecule (1A).  Appropriate iterator
 * is selected based on argument given to setTarget method.  If setTarget method is
 * never called, the all-pair iterator is used by default.<br>
 * This class may be set up to do inter- or intra-species iteration, depending
 * on choice of inner iterators given at construction.
 */
public class ApiMolecule implements AtomsetIteratorMolecule, NearestImageVectorSource {

    /**
     * @param api1A iterator for all pairs formed with a target molecule
     * @param apiAA iterator for all pairs in the phase
     */
    public ApiMolecule(AtomsetIteratorMolecule api1A,
                          AtomsetIteratorPhaseDependent apiAA) {
		this.api1A = api1A;
        this.apiAA = apiAA;
        iterator = apiAA;
	}
    
    public void setTarget(AtomSet targetAtoms) {
        if(targetAtoms == null) {
            iterator = apiAA;
        } else {
            iterator = api1A;
            api1A.setTarget(targetAtoms);
        }
        iterator.setPhase(phase);
    }
    
    public void setPhase(Phase phase) {
        this.phase = phase;
        iterator.setPhase(phase);
    }
       
    public void setDirection(Direction direction) {
        api1A.setDirection(direction);
    }
	/* (non-Javadoc)
	 * @see etomica.AtomIterator#contains(etomica.Atom)
	 */
	public boolean contains(AtomSet atom) {
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
	public AtomSet next() {
		return iterator.next();
	}
	
	public AtomSet peek() {
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
    
    /**
     * Returns the iterator currently being used for 
     * the pair iteration.  If a target had been specified, this will be the
     * Api1A iterator, otherwise it will be the ApiAA iterator.
     */
    public AtomsetIteratorPhaseDependent getCurrentIterator() {
        return iterator;
    }
    
    public Vector getNearestImageVector() {
        return ((NearestImageVectorSource)iterator).getNearestImageVector();
    }
    
    private AtomsetIteratorPhaseDependent iterator;
    private final AtomsetIteratorMolecule api1A; 
    private final AtomsetIteratorPhaseDependent apiAA;
    private Phase phase;

    public AtomsetIteratorMolecule getApi1A() {
        return api1A;
    }
    public AtomsetIteratorPhaseDependent getApiAA() {
        return apiAA;
    }
}
