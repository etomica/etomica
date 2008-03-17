package etomica.modules.pistoncylinder;

import etomica.api.ISimulation;
import etomica.api.IRandom;

import etomica.potential.P1HardMovingBoundary;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;


/**
 * Boundary class for PistonCylinder that accounts for the piston and collision
 * radius when calculating the (available) volume
 * @author andrew
 */
public class BoundaryPistonCylinder extends BoundaryRectangularNonperiodic {

    public BoundaryPistonCylinder(Space _space, IRandom random) {
        super(_space, random);
    }
    
    public void setPistonPotential(P1HardMovingBoundary newPistonPotential) {
        pistonPotential = newPistonPotential;
    }

    public double volume() {
        double collisionDiameter = pistonPotential.getCollisionRadius()*2;
        double v = 1;
        for (int i=0; i<space.D(); i++) {
            if (i == 1) {
                if (space.D() == 2) {
                    // bottom of the box is +dimensions/2, top is wall position
                    v *= (0.5*dimensions.x(i) - pistonPotential.getWallPosition() - collisionDiameter);
                }
                else {
                    // bottom of the box is -dimensions/2, top is wall position
                    v *= (0.5*dimensions.x(i) + pistonPotential.getWallPosition() - collisionDiameter);
                }
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
