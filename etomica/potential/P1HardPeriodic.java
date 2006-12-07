package etomica.potential;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * pseudo-potential for a "collision" time to update colliders for periodic boundaries
 */
 
public class P1HardPeriodic extends Potential1 implements PotentialHard {
    
    public P1HardPeriodic(Space space) {
        super(space);
    }
    
    /**
     * Returns zero.
     */
    public double energy(AtomSet a) {
        return 0.0;
    }
     
    /**
     * Returns zero.
     */
    public double energyChange() {
        return 0.0;
    }
    
    public double collisionTime(AtomSet a, double falseTime) {
        if(!(((Atom)a).getType() instanceof AtomTypeSphere)) {return Double.POSITIVE_INFINITY;}
        Vector v = ((ICoordinateKinetic)((AtomLeaf)a).coord).velocity();
        Vector dim = boundary.getDimensions();
        double tmin = Double.POSITIVE_INFINITY;
        double d2 = 2.0*((AtomTypeSphere)((Atom)a).getType()).getDiameter();
        int D = dim.D();
        for(int i=0; i<D; i++) {
            double t = (dim.x(i)-d2)/v.x(i);
            t = (t < 0) ? -t : t;//abs
            tmin = (t < tmin) ? t : tmin;
        }
        return 0.25*tmin + falseTime;
    }
    
    /**
     * Performs no action.
     */
    public void bump(AtomSet a, double falseTime) { }
    
    /**
     * Returns zero.
     */
    public double lastCollisionVirial() {return 0;}
    
    /**
     * Returns null.
     */
    public Tensor lastCollisionVirialTensor() {return null;}
    
    private static final long serialVersionUID = 1L;
}
   
