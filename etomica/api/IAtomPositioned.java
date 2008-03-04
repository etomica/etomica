package etomica.api;


/**
 * Interface for a IAtom that includes an IVector that defines the atom's
 * position.
 */
public interface IAtomPositioned extends IAtom{
    
    /**
     * Returns the position of the IAtom.  Modifying the returned IVector will
     * alter the IAtom's position.
     */
    public IVector getPosition();

}