package etomica.atom.iterator;

import etomica.api.IAtomLeaf;

/**
 * Interface for an atom iterator that can be altered by setting
 * of an atom.  A neighbor iterator, for example. 
 */
public interface AtomIteratorAtomDependent extends AtomIterator {

	public void setAtom(IAtomLeaf atom);
}
