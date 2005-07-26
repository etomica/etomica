package etomica.space;

import java.io.Serializable;

import etomica.Simulation;
import etomica.Space;

/**
 * Constructs coordinates for spherical atoms, having no orientation dependence.
 * The Simulation instance given at construction provides a Space instance,
 * which determines the spatial dimension of the coordinates (it makes the
 * Vectors that form the coordinates); the Simulation's isDynamic method
 * specifies if the factory returns kinetic or non-kinetic coordinates (kinetic
 * coordinates can hold atom velocities and are needed for molecular dynamics
 * and similar simulations).
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 13, 2005 by kofke
 */
public class CoordinateFactorySphere implements CoordinateFactory, Serializable {

    /**
     * Makes factory using the space and isDynamic field of the given
     * Simulation.
     */
    public CoordinateFactorySphere(Simulation sim) {
        this(sim.space, sim.isDynamic());
    }

    /**
     * @param space
     *            required by Coordinate constructors, and builds the vectors
     *            used by the coordinates to hold positions and velocities.
     * @param isKinetic
     *            flag indicating whether coordinates should include fields for
     *            velocities in addition to positions
     */
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
