package etomica.space;

/**
 * Interface for an atom coordinate that holds vectors for
 * position and velocity.
 */
public interface ICoordinateKinetic extends ICoordinate {

    public IVector getVelocity();
 
}
