/*
 * History
 * Created on Aug 31, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface AtomsetIteratorTargetDependent extends AtomsetIterator {
	public void setTarget(Atom[] targetAtoms);

}