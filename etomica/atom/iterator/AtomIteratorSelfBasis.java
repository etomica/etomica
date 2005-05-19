/*
 * History
 * Created on Sep 8, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;

/**
 * Singleton iterator that takes a basis and returns the basis atom
 * itself as the only iterate.  Targetable so that no atom is given if
 * a target is specified and it is different from the basis.  Expected
 * use is to fill in potential hierarchy for interactions between species
 * having different number of levels forming the molecules (e.g., a single-
 * atom one-layer molecule interacting with a multiatomic).
 */
public class AtomIteratorSelfBasis extends AtomIteratorAdapter implements
		AtomsetIteratorBasisDependent {

	/**
	 * Constructs new iterator in an unset condition.  Must call reset
	 * before beginning iteration.
	 */
	public AtomIteratorSelfBasis() {
		super(new AtomIteratorSinglet());
		singletIterator = (AtomIteratorSinglet)iterator;
	}

	/**
	 * Sets the basis atom, such that the iterator will return the first
	 * atom of the given array.  If array is null or no atom is otherwise
	 * specified, no iterates will be given until a basis atom is specified
	 * via a subsequent call to this method.
	 */
	public void setBasis(AtomSet atoms) {
	    basisAtom = (Atom)atoms;
	}
	
	/**
	 * Puts the iterator in a state so that it is ready to return its iterate.
	 */
	public void reset() {
		if(basisAtom == null || (targetAtom != null && basisAtom != targetAtom)) {
			singletIterator.unset();
		} else {
			singletIterator.setAtom(basisAtom);
			singletIterator.reset();
		}
	}

	/**
	 * Returns 1, indicating that a single-atom basis is expected.
	 */
	public int basisSize() {
		return 1;
	}
    
    public boolean haveTarget(AtomSet target) {
        if(basisAtom == null) {
            return false;
        }
        if(target.count() == 1) {
            return target.getAtom(0) == basisAtom;
        }
        return target.count() == 0;
    }

	/**
	 * Specifies a target atom, such that if the target does not equal the
	 * basis atom at time of reset, no iterates are given.  If a null target
	 * is given, then iteration will always give the most recently specified
	 * basis atom.
	 */
	public void setTarget(AtomSet targetAtoms) {
        switch(targetAtoms.count()) {
        case 0: 
            targetAtom = null;
            break;
        case 1:
            targetAtom = targetAtoms.getAtom(0);
            break;
        default:
            throw new IllegalArgumentException("Can specify at most one target atom to AtomIteratorBasis");
        }
	}

	private Atom basisAtom;
	private Atom targetAtom;
	private final AtomIteratorSinglet singletIterator;
}
