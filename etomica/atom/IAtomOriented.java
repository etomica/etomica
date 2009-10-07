package etomica.atom;

import etomica.api.IAtom;
import etomica.space.IOrientation;

/**
 * Interface for a IAtom that includes an IVector that defines the atom's
 * orientation.
 */
public interface IAtomOriented extends IAtom {
    
    /**
     * Returns the orientation of the IAtom.  Modifying the returned IVector will
     * alter the IAtom's orientation.
     */
    public IOrientation getOrientation();

}