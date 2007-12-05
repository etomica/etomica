package etomica.modules.pistoncylinder;

import etomica.math.geometry.Plane;
import etomica.potential.P1HardMovingBoundary;
import etomica.space.IVector;

/**
 * Wrap a P1HardMovingBoundary and make it look like a Plane.  A boatload of
 * plane methods aren't overriden (there are a lot of them!) and calling them
 * will return garbage (or perhaps even crash).
 * DisplayBoxCanvasG3DSys only calls distanceTo and getD.
 *
 * @author Andrew Schultz
 */
public class PistonPlane extends Plane {
    public PistonPlane(P1HardMovingBoundary pistonPotential) {
        super(pistonPotential.getSpace());
        this.pistonPotential = pistonPotential;
    }
    
    // DisplayBoxCanvasG3DSys calls this
    public double distanceTo(IVector v) {
        return v.x(1) - pistonPotential.getWallPosition();
    }
    
    public double getA() {
        return 0;
    }
    
    public double getB() {
        return 1;
    }
    
    public double getC() {
        return 0;
    }

    // DisplayBoxCanvasG3DSys calls this
    public double getD() {
        return -pistonPotential.getWallPosition();
    }
    
    private static final long serialVersionUID = 1L;
    protected final P1HardMovingBoundary pistonPotential;
}
