/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.ApiInnerVariable;
import etomica.Atom;
import etomica.AtomIteratorBasis;
import etomica.AtomsetIteratorAdapter;
import etomica.AtomsetIteratorBasisDependent;
import etomica.Species;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;

/**
 * Returns iterates from from the childList of a single basis atom.  Behavior is set
 * via iterator directive: if no atom is specified there, all pairs formed from the childList
 * are given; otherwise, if an atom is specified there, pairs will be formed from the
 * childList atoms with the basis' child from which the directive atom is descended.  
 */
public final class ApiIntragroupCell extends AtomsetIteratorAdapter implements
        AtomsetIteratorBasisDependent {

	/**
	 * Constructor makes iterator that must have basis specified and then be 
	 * reset() before iteration.
	 */
	public ApiIntragroupCell(Species[] species) {
		this(new ApiInnerVariable(
                new AtomIteratorBasis(),
                new AtomIteratorNbrCell(species[0], true)));
    }
    
    public ApiIntragroupCell(ApiInnerVariable pairIterator) {   
        super(pairIterator);
		aiOuter = (AtomIteratorBasis)pairIterator.getOuterIterator();
	}

	public void setTarget(Atom[] targetAtoms) {
		aiOuter.setTarget(targetAtoms);
	}

	/**
	 * Specifies the parent atom of the iterates. Pairs given by the
	 * iterator will be formed from the childList of the first atom
	 * in the given array. If argument is null or otherwise does not
	 * specify a childList, iterator will not return iterates 
	 * (hasNext false) until a valid basis is specified.  Length of given
	 * array should match the value returned by setBasis, but if it
	 * is greater no error results; only first atom in array is used.
	 */
	public void setBasis(Atom[] atoms) {
		aiOuter.setBasis(atoms);
	}
	
	/**
	 * Returns the number of atom pairs the iterator will return if
	 * reset and iterated in its present state.
	 */
	public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
	}
	
	/**
	 * Indicates whether the given atom pair will be among the iterates
	 * given by the iterator if reset in its present state.  True only
	 * if an iterated pair would match the atoms as ordered in the given
	 * array.
	 */
	public boolean contains(Atom[] atoms) {
        if(atoms==null || atoms[0]==null || atoms[1]==null || atoms[0]==atoms[1]) return false;
        AtomsetDetect detector = new AtomsetDetect(atoms);
        allAtoms(detector);
        return detector.detectedAtom();
	}

	/**
	 * Returns 1, indicating the the array given to the setBasis method
	 * should have only one element.
	 */
	public int basisSize() {
		return 1;
	}
	
	private final AtomIteratorBasis aiOuter;//local, specifically typed copy
}
