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
public class BoundaryPistonCylinderInner extends BoundaryRectangularNonperiodic {

    public BoundaryPistonCylinderInner(Simulation sim) {
        this(sim.space);
    }

    public BoundaryPistonCylinderInner(Space space) {
        super(space);
    }
    
    public void setPistonPotential(P1HardMovingBoundary newPistonPotential) {
        pistonPotential = newPistonPotential;
    }

    public void setPadding(double newPadding) {
        padding = newPadding;
    }
    
    public double getPadding() {
        return padding;
    }
    
    public double volume() {
        // bottom of the phase is +dimensions/2, top is wall position.
        // subtract padding from each side in each dimension
        return (dimensions.x(0)-2*padding) * (0.5 * dimensions.x(1)-pistonPotential.getWallPosition()-2*padding);
    }

    private P1HardMovingBoundary pistonPotential;
    private double padding;
}
