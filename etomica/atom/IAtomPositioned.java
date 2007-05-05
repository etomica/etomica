package etomica.atom;

import etomica.space.IVector;

/**
 * Interface for a Coordinate that includes a Vector that defines the atom's position.
 */
public interface IAtomPositioned extends IAtom {
    
    /**
     * Returns the position of the IAtom.
     */
    public IVector getPosition();

}