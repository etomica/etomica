package etomica;

/**
 * @author David Kofke
 *
 * AtomPairIterator that wraps an iterator appropriate for all pairs formed from
 * a given atom ("single-loop iteration"), and another for all pairs that can be
 * formed from a given basis ("double-loop iteration").
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
 
public class Api1A extends AtomPairIterator {

	public Api1A(AtomPairIterator api1, AtomPairIterator apiA) {
		this.api1 = api1;
		this.apiA = apiA;
	}
	/**
	 * @see etomica.AtomPairIterator#setBasis(etomica.Atom, etomica.Atom)
	 */
	public void setBasis(Atom a1, Atom a2) {throw new etomica.exception.MethodNotImplementedException();}

	/**
	 * @see etomica.AtomPairIterator#size()
	 */
	public int size() {
		throw new etomica.exception.MethodNotImplementedException();
	}

	/**
	 * @see etomica.AtomPairIterator#hasNext()
	 */
	public boolean hasNext() {
		throw new etomica.exception.MethodNotImplementedException();
	}

	/**
	 * @see etomica.AtomPairIterator#reset(etomica.IteratorDirective)
	 */
	public void reset(IteratorDirective id) {throw new etomica.exception.MethodNotImplementedException();}

	/**
	 * @see etomica.AtomPairIterator#reset()
	 */
	public void reset() {throw new etomica.exception.MethodNotImplementedException();}

	/**
	 * @see etomica.AtomPairIterator#next()
	 */
	public AtomPair next() {
		throw new etomica.exception.MethodNotImplementedException();
	}

	/**
	 * @see etomica.AtomPairIterator#allPairs(etomica.AtomPairAction)
	 */
	public void allPairs(AtomPairAction act) {throw new etomica.exception.MethodNotImplementedException();}

	public void all(Atom basis, IteratorDirective id, AtomPairActive action) {
		if(id.atomCount()==0) apiA.all(basis, id, action);
		else api1.all(basis, id.setSkipFirst(true), action);
	}
	public void all(AtomPair basis, IteratorDirective id, AtomPairActive action) {
		if(id.atomCount()==0) apiA.all(basis, id, action);
		else api1.all(basis, id.setSkipFirst(false), action);
	}
	
	private AtomPairIterator api1 = AtomPairIterator.NULL;
	private AtomPairIterator apiA = AtomPairIterator.NULL;

}
