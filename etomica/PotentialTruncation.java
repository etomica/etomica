package etomica;

/**
 * Class that defines whether and how the interatomic potential is truncated.
 *
 * @author David Kofke
 */
 
public abstract class PotentialTruncation {
    
    protected Potential2Soft parentPotential;
    
    public PotentialTruncation(Potential2Soft potential) {
        parentPotential = potential;
    }
    
    /**
     * Returns true is the truncation makes the potential zero.
     */
    public abstract boolean isZero(double r2);

    /**
     * Transforms the given value of the energy for the truncation applied at
     * the given (square) distance.
     */
    public abstract double uTransform(double r2, double untruncatedValue);

    /**
     * Transforms the given value of the separation-distance derivative for the 
     * truncation applied at the given (square) distance.
     */
    public abstract double duTransform(double r2, double untruncatedValue);

    /**
     * Transforms the given value of the separation-distance second derivative for the 
     * truncation applied at the given (square) distance.
     */
    public abstract double d2uTransform(double r2, double untruncatedValue);
    

    ///************** end of methods for PotentialTruncation ***************
    
    public static final PotentialTruncation NULL = new Null();
    /**
     * No-op version of PotentialTruncation that performs no truncation at all.
     */
     private static final class Null extends PotentialTruncation {
        public Null() {super(null);}
        public boolean isZero(double r2) {return false;}
        public double uTransform(double r2, double untruncatedValue) {return untruncatedValue;}
        public double duTransform(double r2, double untruncatedValue) {return untruncatedValue;}
        public double d2uTransform(double r2, double untruncatedValue) {return untruncatedValue;}
     }//end of Null
}//end of PotentialTruncation