/*
 * Created on Aug 5, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface AtomPairPhaseIterator extends AtomPairIterator {
	public void all(Phase phase, IteratorDirective id, AtomPairActive action);
	
	public void setBasis(Phase phase);
	
}
