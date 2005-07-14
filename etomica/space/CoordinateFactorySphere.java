package etomica.space;

import etomica.Space;

/**
 * Constructs coordinates for spherical atoms, having no orientation dependence.
 * Can be configured to return kinetic or non-kinetic coordinates (kinetic
 * coordinates can hold atom velocities and are needed for molecular dynamics
 * and similar simulations).  The Space instance given at construction determines
 * the spatial dimension of the coordinates; it makes the Vectors that form the 
 * coordinates.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 13, 2005 by kofke
 */
public class CoordinateFactorySphere implements CoordinateFactory {

    /**
     * Makes factory with default value of isKinetic = true.
     * 
     * @param space
     *            required by Coordinate constructors, and builds the vectors
     *            used by the coordinates to hold positions and velocities.
     */
    public CoordinateFactorySphere(Space space) {
        this(space, true);
    }

    public CoordinateFactorySphere(Space space, boolean isKinetic) {
        this.space = space;
        this.isKinetic = isKinetic;
    }

    /**
     * Returns a new instance of Coordinate or CoordinateKinetic, according to
     * the current value of isKinetic.
     */
    public ICoordinate makeCoordinate() {
        return isKinetic ? new CoordinateKinetic(space) : new Coordinate(space);
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
