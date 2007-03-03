package etomica.space;

/**
 * Interface for an atom coordinate that holds vectors for
 * position and velocity.
 */

/*
 * History
 * Created on Jan 25, 2005 by kofke
 */
public interface ICoordinateKinetic extends ICoordinate {

    public IVectorRandom getVelocity();
 
}
