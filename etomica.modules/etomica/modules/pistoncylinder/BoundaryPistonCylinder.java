package etomica.modules.pistoncylinder;

import etomica.potential.P1HardMovingBoundary;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.util.IRandom;

/**
 * Boundary class for PistonCylinder that accounts for the piston and collision
 * radius when calculating the (available) volume
 * @author andrew
 */
public class BoundaryPistonCylinder extends BoundaryRectangularNonperiodic {


    public BoundaryPistonCylinder(Simulation sim) {
        this(sim.getSpace(), sim.getRandom());
    }

    public BoundaryPistonCylinder(Space space, IRandom random) {
        super(space, random);
    }
    
    public void setPistonPotential(P1HardMovingBoundary newPistonPotential) {
        pistonPotential = newPistonPotential;
    }

    public double volume() {
        double collisionDiameter = pistonPotential.getCollisionRadius()*2;
        // bottom of the phase is +dimensions/2, top is wall position
        return (dimensions.x(0) - collisionDiameter) * 
               (0.5*dimensions.x(1) - pistonPotential.getWallPosition() - collisionDiameter);
    }

    private P1HardMovingBoundary pistonPotential;
    private static final long serialVersionUID = 1L;
}
