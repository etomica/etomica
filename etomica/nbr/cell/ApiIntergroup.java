/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;
import etomica.*;
import etomica.lattice.BravaisLattice;
import etomica.math.geometry.Polyhedron;
import etomica.nbr.cell.IteratorFactoryCell.CellSequencer;
import etomica.utility.java2.HashMap;

/**
 * Iterator that returns pairs formed using two different basis atoms, so 
 * that the iterates are taken from two atom different groups.  
 */
public final class ApiIntergroup extends AtomsetIteratorAdapter implements
		AtomsetIteratorBasisDependent {

	public ApiIntergroup() {
		super(new ApiInnerFixed(
				new AtomIteratorBasis(),
				new AtomIteratorBasis()));
		ApiInnerFixed pairIterator = (ApiInnerFixed)iterator;
		aiOuter = (AtomIteratorBasis)pairIterator.getOuterIterator();
		aiInner = (AtomIteratorBasis)pairIterator.getInnerIterator();
	}
	
	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorBasisDependent#setDirective(etomica.IteratorDirective)
	 */
	public void setTarget(Atom[] targetAtoms) {
		aiOuter.setTarget(targetAtoms);
	}

	/**
	 * Specifies the basis, which identifies the atoms subject
	 * to iteration. The given array should be of length 2 (at least);
	 * first atom in array specifies the basis for the outer-loop
	 * iteration, and second atom specifies the basis for the inner-loop
	 * iteration.  In each case, if the basis atom is not a leaf atom, its children will
	 * be the subject of iteration.  If the basis atom is a leaf, it
	 * will itself be the iterate.  If given array is null, or if its length
	 * is not at least 2, iterator will give no iterates until a proper basis
	 * is specified via another call to this method.
	 */
	public void setBasis(Atom[] basisAtoms) {
		if(basisAtoms == null || basisAtoms.length < 2) {
			aiOuter.setBasis(null);
		} else {
			atom[0] = basisAtoms[0];
			aiOuter.setBasis(atom);
			atom[0] = basisAtoms[1];
			aiInner.setBasis(atom);
		}
	}

	/**
	 * Returns 2, indicating that the setBasis method expects an 
	 * array of two atoms.
	 */
	public int basisSize() {
		return 2;
	}

	private final AtomIteratorBasis aiOuter;
	private final AtomIteratorBasis aiInner;
	private final Atom[] atom = new Atom[1];

}
