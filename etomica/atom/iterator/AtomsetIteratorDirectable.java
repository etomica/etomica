/*
 * History
 * Created on Sep 3, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.AtomsetIterator;
import etomica.IteratorDirective;

/**
 * 
 * Interface for an iterator that can interpret specification of
 * direction UP or DOWN.
 *
 * @author David Kofke
 *
 */
public interface AtomsetIteratorDirectable extends AtomsetIterator {

	public void setDirection(IteratorDirective.Direction direction);
	
}
