package etomica;

/**
 * Class that defines whether and how the interatomic potential is truncated.
 *
 * @author David Kofke
 */
 
public abstract class PotentialTruncation {
    
    protected final Potential2 potential;
    
    public PotentialTruncation(Potential2 potential) {
        this.potential = potential;
        if(potential != null) {//expect null only for inner Null class defined below
            Potential0GroupLrc lrcMaster = potential.parentSimulation().hamiltonian.potential.lrcMaster();
            makeLrcPotential(lrcMaster); //adds this to lrcMaster
        }
    }
    
    /**
     * Returns true if the truncation makes the potential zero at the given separation.
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
    
    /**
     * Returns a class that calculates the long-range contribution to the potential
     * that becomes neglected by the truncation.  Assumes a uniform distribution
     * of atoms beyond this truncation's cutoff distance.
     */
    public abstract Potential0Lrc makeLrcPotential(PotentialGroup parent);
    

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
        public Potential0Lrc makeLrcPotential(PotentialGroup parent) {return null;}
     }//end of Null
}//end of PotentialTruncation