package etomica.modules.pistoncylinder;

import etomica.potential.P1HardMovingBoundary;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;

/**
 * Boundary class for PistonCylinder that accounts for the piston when 
 * calculating the volume
 * @author andrew
 */
public class BoundaryPistonCylinder extends BoundaryRectangularNonperiodic {

    public BoundaryPistonCylinder(Simulation sim) {
        this(sim.space);
    }

    public BoundaryPistonCylinder(Space space) {
        super(space);
    }
    
    public void setPistonPotential(P1HardMovingBoundary newPistonPotential) {
        pistonPotential = newPistonPotential;
    }

    public double volume() {
        // bottom of the phase is +dimensions/2, top is wall position
        return dimensions.x(0) * (0.5 * dimensions.x(1)-pistonPotential.getWallPosition());
    }

    private P1HardMovingBoundary pistonPotential;
}
