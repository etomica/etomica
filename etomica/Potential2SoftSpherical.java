package etomica;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 * Truncation (if defined for the instance) is applied in the energy, virial, 
 * hypervirial and gradient methods.  Subclasses must provide concrete definitions
 * for the energy (method u(double)) and its derivatives without applying truncation.
 *
 * @author David Kofke
 */
 
 /* History
  * (this change was not implemented) 10/12/02 (DAK) modified virial method to reflect change in definition of virial in Potential2Soft
  * 08/11/03 modified gradient method to colsolidate two vector operations
  * into one method call
  */
public abstract class Potential2SoftSpherical extends Potential2 implements Potential2.Soft {
   
   private final Space.Vector work1;
   private final double rD;// = 1/D
   
   public Potential2SoftSpherical(SimulationElement parent) {
        super(parent);
        rD = 1.0/(double)simulation().space.D();
        work1 = simulation().space().makeVector();
   }
   public Potential2SoftSpherical(SimulationElement parent, PotentialTruncation trunc) {
     //constructors repeat code rather than call the other because superclass constructors
     //define truncation differently.  Since truncation field is final it cannot be
     //subsequently changed
        super(parent, trunc);
        rD = 1.0/(double)simulation().space.D();
        work1 = simulation().space.makeVector();
   }
        
   
   /**
    * The pair energy u(r^2) with no truncation applied.
    * @param the square of the distance between the particles.
    */
    public abstract double u(double r2);
        
   /**
    * The derivative of the pair energy, times the separation r: r du/dr.
    * No truncation applied.
    */
    public abstract double du(double r2);
        
   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.  No truncation applied.
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
    
    /**
     * Energy of the pair as given by the u(double) method, with application
     * of any PotentialTruncation that may be defined for the potential.
     */
    public double energy(Atom[] pair) {
     //   return u(pair.r2());
    	cPair.reset(pair[0].coord,pair[1].coord);
    	double r2 = cPair.r2();
        if(potentialTruncation.isZero(r2)) return 0.0;
        else return potentialTruncation.uTransform(r2, u(r2));
    }
    
    /**
     * Virial of the pair as given by the du(double) method, with application
     * of any PotentialTruncation that may be defined for the potential.
     */
    public double virial(Atom[] pair) {
    	cPair.reset(pair[0].coord,pair[1].coord);
        double r2 = cPair.r2();
        if(potentialTruncation.isZero(r2)) return 0.0;
        else return potentialTruncation.duTransform(r2, du(r2));
    }
    
    /**
     * Hypervirial of the pair as given by the du(double) and d2u(double) methods, with application
     * of any PotentialTruncation that may be defined for the potential.
     */
    public double hyperVirial(Atom[] pair) {
    	cPair.reset(pair[0].coord,pair[1].coord);
        double r2 = cPair.r2();
        if(potentialTruncation.isZero(r2)) return 0.0;
        else return potentialTruncation.d2uTransform(r2, d2u(r2)) + potentialTruncation.duTransform(r2, du(r2));
    }
    
    /**
     * Gradient of the pair potential as given by the du(double) method, with application
     * of any PotentialTruncation that may be defined for the potential.
     */
    public Space.Vector gradient(Atom[] pair) {
  //  	System.out.println(((P2LennardJones)this).getSigma()+"  "+((AtomType.Sphere)pair.atom1().type).diameter);
    	cPair.reset(pair[0].coord,pair[1].coord);
        double r2 = cPair.r2();
        if(potentialTruncation.isZero(r2)) work1.E(0.0);
        else {
            double v = potentialTruncation.duTransform(r2, du(r2));
            work1.Ea1Tv1(v/r2,cPair.dr());
//            work1.E(pair.dr());
//            work1.TE(v/r2);
        }
        return work1;
    }
    
    /**
     * Same as uInt.
     */
    public double integral(double rC) {
        return uInt(rC);
    }
}//end of Potential2SoftSpherical
