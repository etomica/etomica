package etomica.potential;

import etomica.Atom;
import etomica.AtomSet;
import etomica.Space;
import etomica.atom.AtomTypeSphere;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * pseudo-potential for a "collision" time to update colliders for periodic boundaries
 */
 
public class P1HardPeriodic extends Potential1 implements PotentialHard {
    
    public P1HardPeriodic(Space space) {
        super(space);
    }
    
    public double energy(AtomSet a) {
        return 0.0;
    }
     
    public double energyChange() {
        return 0.0;
    }
    
    public double collisionTime(AtomSet a, double falseTime) {
        if(!(((Atom)a).type instanceof AtomTypeSphere)) {return Double.MAX_VALUE;}
        Vector v = ((ICoordinateKinetic)((Atom)a).coord).velocity();
        Vector dim = boundary.dimensions();
        double tmin = Double.MAX_VALUE;
        double d2 = 2.0*((AtomTypeSphere)((Atom)a).type).diameter((Atom)a);
        int D = dim.D();
        for(int i=0; i<D; i++) {
            double t = (dim.x(i)-d2)/v.x(i);
            t = (t < 0) ? -t : t;//abs
            tmin = (t < tmin) ? t : tmin;
        }
        return 0.25*tmin + falseTime;
    }
                
    public void bump(AtomSet a, double falseTime) { }
        
    public double lastCollisionVirial() {return 0;}
    
    public Tensor lastCollisionVirialTensor() {return null;}
    
}
   
