package etomica.atom.iterator;

import etomica.api.IMolecule;

/**
 * Interface for an atom iterator that can be altered by setting
 * of an atom.  A neighbor iterator, for example. 
 */
public interface MoleculeIteratorMoleculeDependent extends MoleculeIterator {

	public void setMolecule(IMolecule atom);
}
