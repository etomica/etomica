package etomica.space;

/**
 * Interface for a dynamical coordinate that includes an orientation in
 * addition to a spatial position.  
 */
/*
 * History
 * Created on Jan 25, 2005 by kofke
 */
public interface ICoordinateAngularKinetic extends ICoordinateKinetic, ICoordinateAngular {
    
    public Vector angularVelocity(); //angular velocity vector in space-fixed frame
    
}