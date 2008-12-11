package etomica.atom.iterator;

import etomica.api.IBox;

/**
 * Interface indicating that an iterator can determine appropriate
 * atoms for iteration given an arbitrary box.  
 */
public interface MoleculesetIteratorBoxDependent extends MoleculesetIterator {

    /**
     * Sets the Box to pull iterates from
     * @throws a NullPointerException if the Box is null
     */
	public void setBox(IBox box);
	
}
