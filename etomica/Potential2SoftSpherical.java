package etomica;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 */
public abstract class Potential2SoftSpherical extends Potential2 implements Potential2Soft {
   
   public Potential2SoftSpherical(Simulation sim) {
        super(sim);
   }
   
   /**
    * The pair energy u(r^2).
    * @param the square of the distance between the particles.
    */
    public abstract double u(double r2);
        
   /**
    * The derivative of the pair energy, times the separation r: r du/dr.
    */
    public abstract double du(double r2);
        
   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public abstract double d2u(double r2);
        
   /**
    * Long-range contribution to the energy, calculated assuming a uniform distribution
    * of other particles from the cutoff to infinity.
    */
    public abstract double uLRC();
        
   /**
    * Long-range contribution to the energy first derivative.
    */
    public abstract double duLRC();
        
   /**
    * Long-range contribution to the energy second derivative.
    */
    public abstract double d2uLRC();
    
    public double energy(AtomPair pair) {
        return u(pair.r2());
    }
    
    public double virial(AtomPair pair) {
        return du(pair.r2());
    }
    
    public double hyperVirial(AtomPair pair) {
        double r2 = pair.r2();
        return d2u(r2) + du(r2);
    }
    
    public Space.Vector gradient(AtomPair pair) {
        double r2 = pair.r2();
        double v = du(r2);
        work1.E(pair.dr());
        work1.TE(v/r2);
        return work1;
    }
    
}//end of Potential2SoftSpherical
