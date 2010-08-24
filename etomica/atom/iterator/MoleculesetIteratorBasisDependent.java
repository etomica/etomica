package etomica.atom.iterator;

import etomica.api.IMolecule;
import etomica.api.IMoleculeList;

/**
 * Interface for an MoleculeIterator that can be conditioned with
 * target and basis molecules.
 */
public interface MoleculesetIteratorBasisDependent extends MoleculesetIteratorTargetable {

	/**
	 * 
	 */
    public void setBasis(IMoleculeList molecules);
    
    /**
     * Indicates the size of the basis needed to set the iterator.
     * The length of the array given to setBasis should be this value.
     * @return the size of the basis for this iterator.
     */
    public int basisSize();
    
    /**
     * Returns true if the iterator with its current basis 
     * would return an iterate for the given target.
     */
    public boolean haveTarget(IMolecule target);

}
