package etomica.potential;

import etomica.api.IAtomKinetic;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomTypeSphere;
import etomica.api.IVectorMutable;
import etomica.api.IVector;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * pseudo-potential for a "collision" time to update colliders for periodic boundaries
 */
 
public class P1HardPeriodic extends Potential1 implements PotentialHard {

    /**
     * Returns an instance of P1HardPeriodic with sigma = NaN.  call setSigma
     * to set the value you want.
     */
    public P1HardPeriodic(ISpace space) {
        this(space, Double.NaN);
        // use NaN so they'll have to call setSigma later
    }

    /**
     * Returns an instance of P1HardPeriodic with the given value of sigma (the
     * maximum distance between two atoms where they interact)
     */
    public P1HardPeriodic(ISpace space, double sigma) {
        super(space);
        this.sigma = sigma;
    }
    
    /**
     * Returns zero.
     */
    public double energy(IAtomList a) {
        return 0.0;
    }
     
    /**
     * Returns zero.
     */
    public double energyChange() {
        return 0.0;
    }
    
    public double collisionTime(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.getAtom(0);
        if(!(((IAtom)atom).getType() instanceof IAtomTypeSphere)) {return Double.POSITIVE_INFINITY;}
        IVectorMutable v = atom.getVelocity();
        IVector dim = boundary.getBoxSize();
        double tmin = Double.POSITIVE_INFINITY;
        double d2 = 2.0*sigma;
        int D = dim.getD();
        for(int i=0; i<D; i++) {
            double t = (dim.getX(i)-d2)/v.getX(i);
            t = (t < 0) ? -t : t;//abs
            tmin = (t < tmin) ? t : tmin;
        }
        return 0.25*tmin + falseTime;
    }
    
    public void setSigma(double newSigma) {
        sigma = newSigma;
    }
    
    public double getSgima() {
        return sigma;
    }
    
    /**
     * Performs no action.
     */
    public void bump(IAtomList a, double falseTime) { }
    
    /**
     * Returns zero.
     */
    public double lastCollisionVirial() {return 0;}
    
    /**
     * Returns null.
     */
    public Tensor lastCollisionVirialTensor() {return null;}
    
    private static final long serialVersionUID = 1L;
    protected double sigma;
}
   
