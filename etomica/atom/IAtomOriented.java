package etomica.atom;

import etomica.api.IAtomPositioned;
import etomica.space.IOrientation;

/**
 * Interface for a IAtom that includes an IVector that defines the atom's
 * position.
 */
public interface IAtomOriented extends IAtomPositioned {
    
    /**
     * Returns the orientation of the IAtom.  Modifying the returned IVector will
     * alter the IAtom's position.
     */
    public IOrientation getOrientation();

}