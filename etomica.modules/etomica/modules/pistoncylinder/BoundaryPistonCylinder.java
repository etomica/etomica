package etomica.modules.pistoncylinder;

import etomica.potential.P1HardMovingBoundary;
import etomica.simulation.ISimulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.util.IRandom;

/**
 * Boundary class for PistonCylinder that accounts for the piston and collision
 * radius when calculating the (available) volume
 * @author andrew
 */
public class BoundaryPistonCylinder extends BoundaryRectangularNonperiodic {


    public BoundaryPistonCylinder(ISimulation sim) {
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
        double v = 1;
        for (int i=0; i<space.D(); i++) {
            if (i == 1) {
                // bottom of the box is +dimensions/2, top is wall position
                v *= (0.5*dimensions.x(i) - pistonPotential.getWallPosition() - collisionDiameter);
            }
            else {
                v *= dimensions.x(i) - collisionDiameter;
            }
        }
        return v;
    }

    private P1HardMovingBoundary pistonPotential;
    private static final long serialVersionUID = 1L;
}
