package etomica.modules.osmosis;
import etomica.EtomicaInfo;
import etomica.atom.AtomSet;
import etomica.atom.IAtomKinetic;
import etomica.potential.Potential1;
import etomica.potential.PotentialHard;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.Length;

/**
 */
 
public class P1HardWall extends Potential1 implements PotentialHard {
    
    private static final long serialVersionUID = 1L;
    private double collisionRadius;
    
    public P1HardWall(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().atomSize);
    }
    
    public P1HardWall(Space space, double sigma) {
        super(space);
        collisionRadius = sigma;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Harmonic potential at the phase boundaries");
        return info;
    }

    public double energy(AtomSet a) {
        double e = 0.0;
        //XXX ignore atoms in the wall.  this can happen due to bogus initial configurations
//        if (Math.abs(((AtomLeaf)a).coord.position().x(0)) < collisionRadius) {
//            e = Double.MAX_VALUE;
//        }
        return e;
    }

     
    public double collisionTime(AtomSet a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a;
        IVector r = atom.getPosition();
        IVector v = atom.getVelocity();
        double vx = v.x(0);
        // cheat!  We want to ignore collisions from the left.  The initial
        // config sometimes plops atoms on the left.
        if (vx > 0) {
            return Double.POSITIVE_INFINITY;
        }
        double rx = r.x(0) + vx * falseTime;
        double t = (vx > 0.0) ? - collisionRadius : collisionRadius;
        t = (t - rx) / vx;
        if (t < 0) {
            // moving away from the wall
            t = Double.POSITIVE_INFINITY;
        }
        return t+falseTime;
    }

    public void bump(AtomSet a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a;
        IVector v = atom.getVelocity();

        v.setX(0,-v.x(0));

        double newP = atom.getPosition().x(0) - falseTime*v.x(0)*2.0;
        atom.getPosition().setX(0,newP);
    }

    public double energyChange() {
        return 0;
    }
    
    /**
     * not yet implemented
     */
    public double lastCollisionVirial() {return Double.NaN;}
    
    /**
     * not yet implemented.
     */
    public Tensor lastCollisionVirialTensor() {return null;}
    
        
    /**
     * Distance from the center of the sphere to the boundary at collision.
     */
    public void setCollisionRadius(double d) {collisionRadius = d;}
    /**
     * Distance from the center of the sphere to the boundary at collision.
     */
    public double getCollisionRadius() {return collisionRadius;}
    /**
     * Indicates collision radius has dimensions of Length.
     */
    public etomica.units.Dimension getCollisionRadiusDimension() {return Length.DIMENSION;}

}
   
