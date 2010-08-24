package etomica.atom.iterator;


/**
 * 
 * Interface for an iterator that can interpret specification of
 * direction UP or DOWN.
 *
 * @author Tai Boon Tan
 */
public interface MoleculesetIteratorDirectable extends MoleculesetIterator {

	public void setDirection(IteratorDirective.Direction direction);
	
}
