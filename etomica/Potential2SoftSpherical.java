package etomica;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 *
 * @author David Kofke
 */
public abstract class Potential2SoftSpherical extends Potential2Soft {
   
   public static String VERSION = "Potential2SoftSpherical:01.07.05/"+Potential2.VERSION;
   
   public Potential2SoftSpherical(PotentialGroup parent) {
        super(parent);
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
    * Integral of the potential, used to evaluate corrections for potential truncation.
    * Specifically, this is the integral from rC (the argument) to infinity of
    * u(r) A r^(D-1), where D is the spatial dimension, and A is the area of a unit
    * sphere in D dimensions.  Normally, the long-range potential correction would be obtained
    * by multiplying this quantity by the pair density nPairs/V, where nPairs is the number of pairs of atoms
    * affected by this potential, and V is the volume they occupy.
    */
    public abstract double uInt(double rC);
      
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
    
    public double integral(double rC) {
        return uInt(rC);
    }
}//end of Potential2SoftSpherical
