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
 
public final class Api1A implements AtomPairIterator {

	public Api1A(AtomPairIterator api1, AtomPairIterator apiA) {
		this.api1 = api1;
		this.apiA = apiA;
	}
	/**
	 * @see etomica.AtomPairIterator#setBasis(etomica.Atom, etomica.Atom)
	 */
	public void setBasis(Atom a1, Atom a2) {
		basis1 = a1;
		basis2 = a2;
	}

	/**
	 * @see etomica.AtomPairIterator#size()
	 */
	public int size() {
		return api.size();
	}

	/**
	 * @see etomica.AtomPairIterator#hasNext()
	 */
	public boolean hasNext() {
		return api.hasNext();
	}

	/**
	 * @see etomica.AtomPairIterator#reset(etomica.IteratorDirective)
	 */
	public void reset(IteratorDirective id) {
		api = (id.atomCount() == 0) ? apiA : api1;
		api.setBasis(basis1, basis2);
		api.reset(id);
	}

	/**
	 * @see etomica.AtomPairIterator#reset()
	 */
	public void reset() {throw new etomica.exception.MethodNotImplementedException();}

	/**
	 * @see etomica.AtomPairIterator#next()
	 */
	public AtomPair next() {
//		AtomPair nextPair = api.next(); System.out.println(nextPair.atom1.toString()+nextPair.atom2.toString()); return nextPair; 
		return api.next();
	}

	/**
	 * @see etomica.AtomPairIterator#allPairs(etomica.AtomPairAction)
	 */
	public void allPairs(AtomPairAction act) {throw new etomica.exception.MethodNotImplementedException();}

	public void all(AtomSet basis, IteratorDirective id, final AtomSetActive action) {
		if(basis == null || !(action instanceof AtomPairActive)) return;
		switch(basis.nBody()) {
			case 1: all((Atom)basis, id, (AtomPairActive)action); break;
			case 2: all((AtomPair)basis, id, (AtomPairActive)action); break;
		}
	}
	public void all(Atom basis, IteratorDirective id, AtomPairActive action) {
		if(id.atomCount()==0) apiA.all(basis, id, action);
		else api1.all(basis, id.setSkipFirst(true), action);
	}
	public void all(AtomPair basis, IteratorDirective id, AtomPairActive action) {
		if(id.atomCount()==0) apiA.all(basis, id, action);
		else api1.all(basis, id.setSkipFirst(false), action);
	}
	
	private final AtomPairIterator api1;
	private final AtomPairIterator apiA;
	private AtomPairIterator api;
	private Atom basis1, basis2;

}
