package etomica.potential;

import etomica.Atom;
import etomica.Space;
import etomica.atom.AtomTypeSphere;
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space.Boundary.Periodic;

/**
 * pseudo-potential for a "collision" time to update colliders for periodic boundaries
 */
 
public class P1HardPeriodic extends Potential1 implements PotentialHard {
    
    public P1HardPeriodic(Space space) {
        super(space);
    }
    
    public double energy(Atom[] a) {
        return 0.0;
    }
     
    public double energyChange() {
        return 0.0;
    }
    
    public double collisionTime(Atom[] a, double falseTime) {
        if(boundary instanceof Boundary.Periodic) {
            if(!(a[0].type instanceof AtomTypeSphere)) {return Double.MAX_VALUE;}
            Vector p = a[0].coord.momentum();
            Vector dim = boundary.dimensions();
            double tmin = Double.MAX_VALUE;
            double d2 = 2.0*((AtomTypeSphere)a[0].type).diameter(a[0]);
            int D = dim.D();
            for(int i=0; i<D; i++) {
                double t = (dim.x(i)-d2)/p.x(i);
                t = (t < 0) ? -t : t;//abs
                tmin = (t < tmin) ? t : tmin;
            }
            return 0.25*a[0].coord.mass()*tmin + falseTime; //0.5*m*min of (dim.x/p.x, dim.y/p.y, etc.)
      //      return 0.25*atom.mass()*dim.D(p).abs().min(); //0.5*m*min of (dim.x/p.x, dim.y/p.y, etc.)
            //assumes range of potential is .le. diameter, simulation box is square (or x is smaller dimension)
        }
        return Double.MAX_VALUE;
    }
                
    public void bump(Atom[] a, double falseTime) { }
        
    public double lastCollisionVirial() {return 0;}
    
    public Tensor lastCollisionVirialTensor() {return null;}
    
}
   
