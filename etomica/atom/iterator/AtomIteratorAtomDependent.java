/*
 * History
 * Created on Aug 12, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomIterator;

/**
 * Interface for an atom iterator that can be altered by setting
 * of an atom.  A neighbor iterator, for example. 
 */
public interface AtomIteratorAtomDependent extends AtomIterator {

	public void setAtom(Atom atom);
}
