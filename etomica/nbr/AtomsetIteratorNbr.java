/*
 * History
 * Created on Sep 22, 2004 by kofke
 */
package etomica.nbr;

import etomica.AtomsetIterator;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface AtomsetIteratorNbr extends AtomsetIterator {

	/**
	 * 
	 */
	public NeighborCriterion nbrCriterion();
}
