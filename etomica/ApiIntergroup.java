/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ApiIntergroup extends AtomsetIteratorAdapter implements
		AtomsetIteratorDirectable {

	public ApiIntergroup() {
		super(new ApiInnerVariable(
				new AtomIteratorDirectable(),
				new AtomIteratorSequencerList()));
		pairIterator = (ApiInnerFixed)iterator;
		aiOuter = (AtomIteratorDirectable)pairIterator.getOuterIterator();
		aiInner = (AtomIteratorListSimple)pairIterator.getInnerIterator();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#setDirective(etomica.IteratorDirective)
	 */
	public void setDirective(IteratorDirective id) {
		aiOuter.setDirective(id);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#setBasis(etomica.Atom[])
	 */
	public void setBasis(Atom[] atoms) {
		atom[0] = atoms[0];
		aiOuter.setBasis(atom);
		aiInner.setList(((AtomTreeNodeGroup)atoms[1].node).childList);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#basisSize()
	 */
	public int basisSize() {
		return 2;
	}

	private final ApiInnerFixed pairIterator;
	private final AtomIteratorDirectable aiOuter;
	private final AtomIteratorListSimple aiInner;
	private final Atom[] atom = new Atom[1];

}
