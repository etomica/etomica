/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica;

import etomica.exception.MethodNotImplementedException;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ApiIntragroup extends AtomsetIteratorAdapter implements
		AtomsetIteratorDirectable {

	/**
	 * @param iterator
	 */
	public ApiIntragroup() {
		super(new ApiInnerVariable(
					new AtomIteratorDirectable(),
					new AtomIteratorSequencerList()));
		pairIterator = (ApiInnerVariable)iterator;
		aiOuter = (AtomIteratorDirectable)pairIterator.getOuterIterator();
		aiInner = (AtomIteratorSequencerList)pairIterator.getInnerIterator();
		aiInner.setSkippingFirst(true);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#setDirective(etomica.IteratorDirective)
	 */
	public void setDirective(IteratorDirective id) {
		aiOuter.setDirective(id);
		doBoth = id.atomCount() == 1;
	}
	
	public void reset() {
		if(doBoth) {
			upListNow = true;
			aiInner.setIterationDirection(IteratorDirective.UP);
			pairIterator.reset();
			if(!pairIterator.hasNext()) {
				aiInner.setIterationDirection(IteratorDirective.DOWN);
				upListNow = false;
				pairIterator.reset();
			}
		}		
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#setBasis(etomica.Atom[])
	 */
	public void setBasis(Atom[] atoms) {
		basis = atoms[0];
		aiOuter.setBasis(atoms);
	}
	
	public Atom[] next() {
		Atom[] next = pairIterator.next();
		if(upListNow && doBoth && !pairIterator.hasNext()) {
			aiInner.setIterationDirection(IteratorDirective.DOWN);
			pairIterator.reset();
			upListNow = false;
		}
		return next;
	}
	
	public int size() {
		int n = ((AtomTreeNodeGroup)basis.node).childList.size();
		return doBoth ? n-1 : n*(n-1)/2;
	}
	
	//TODO finish this
	public boolean contains(Atom[] atoms) {
//		AtomList list = ((AtomTreeNodeGroup)basis.node).childList;
//		if(doBoth) return (atoms[0] == basis)
//		else return (atoms[0] != atoms[1]) && list.contains(atoms[0]) && list.contains(atoms[1]);
		throw new MethodNotImplementedException("not yet");
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#basisSize()
	 */
	public int basisSize() {
		return 1;
	}
	
	private final ApiInnerVariable pairIterator;
	private final AtomIteratorDirectable aiOuter;
	private final AtomIteratorSequencerList aiInner;
	private boolean doBoth, upListNow;
	private Atom basis;
}
