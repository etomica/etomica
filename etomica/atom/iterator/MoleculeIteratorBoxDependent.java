package etomica.atom.iterator;


/**
 * Interface indicating that an atom iterator can determine appropriate
 * atoms for iteration given an arbitrary box.  
 */
public interface MoleculeIteratorBoxDependent extends MoleculeIterator, MoleculesetIteratorBoxDependent {

}
