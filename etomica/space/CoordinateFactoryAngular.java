package etomica.space;

import etomica.simulation.Simulation;

/**
 * Constructs coordinates for atoms that have an orientation. The Simulation
 * instance given at construction provides a Space instance, which determines
 * the spatial dimension of the coordinates (it makes the Vectors and
 * Orientations that form the coordinates); the Simulation's isDynamic method
 * specifies if the factory returns kinetic or non-kinetic coordinates (kinetic
 * coordinates can hold atom velocities and are needed for molecular dynamics
 * and similar simulations).
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
     * Makes factory using the space and isDynamic field of the given
     * Simulation.
     */
    public CoordinateFactoryAngular(Simulation sim) {
        this(sim.getSpace(), sim.isDynamic());
    }

    /**
     * @param space
     *            required by Coordinate constructors, and builds the vectors
     *            and orientations used by the coordinates to hold positions and
     *            velocities.
     * @param isKinetic
     *            flag indicating whether coordinates should include fields for
     *            linear and angular velocities in addition to positions and
     *            orientations
     */
    public CoordinateFactoryAngular(Space space, boolean isKinetic) {
        this.space = space;
        this.isKinetic = isKinetic;
    }

    /**
     * Returns a new instance of CoordinateAngular or CoordinateAngularKinetic,
     * according to the current value of isKinetic.
     */
    public ICoordinate makeCoordinate() {
        return isKinetic ? new CoordinateAngularKinetic(space)
                : new CoordinateAngular(space);
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
