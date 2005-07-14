package etomica.space;

import etomica.Space;

/**
 * Constructs coordinates for atoms that have an orientation.
 * Can be configured to return kinetic or non-kinetic coordinates (kinetic
 * coordinates can hold atom velocities and are needed for molecular dynamics
 * and similar simulations).  The Space instance given at construction determines
 * the spatial dimension of the coordinates; it makes the Vectors and Orientations 
 * that form the coordinates.
 * 
 * @author David Kofke
 *  
 */


//TODO provide means to set the type of angular coordinate

/*
 * History Created on Jul 13, 2005 by kofke
 */
public class CoordinateFactoryAngular implements CoordinateFactory {

    /**
     * Makes factory with default value of isKinetic = true.
     * 
     * @param space
     *            required by Coordinate constructors, and builds the vectors
     *            used by the coordinates to hold positions and velocities.
     */
    public CoordinateFactoryAngular(Space space) {
        this(space, true);
    }

    public CoordinateFactoryAngular(Space space, boolean isKinetic) {
        this.space = space;
        this.isKinetic = isKinetic;
    }

    /**
     * Returns a new instance of CoordinateAngular or CoordinateAngularKinetic, according to
     * the current value of isKinetic.
     */
    public ICoordinate makeCoordinate() {
        return isKinetic ? new CoordinateAngularKinetic(space) : new CoordinateAngular(space);
    }

    /**
     * Mutator method for isKinetic flag.
     */
    public void setKinetic(boolean kinetic) {
        isKinetic = kinetic;
    }

    /**
     * Accessor method for isKinetic flag.
     */
    public boolean isKinetic() {
        return isKinetic;
    }

    private final Space space;
    private boolean isKinetic;

}
