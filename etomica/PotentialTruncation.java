package etomica;

import etomica.Space.CoordinatePair;

/**
 * Class that defines whether and how the interatomic potential is truncated.
 *
 * @see Potential0Lrc
 * @author David Kofke
 */
 
 /* History of changes
  * 07/13/02 (DAK) Restructured instantiation of LRC potential
  * 01/24/03 (DAK) Added isZero(AtomSet) method as part of Potential redesign
  */
 
public abstract class PotentialTruncation {
    
    public PotentialTruncation() {}
    
    /**
     * Truncation used by potential group.
     * @see PotentialCalculation.PotentialGroupWrapper
     * @param atomSet  Candidate set of atoms for truncation
     * @return boolean true if potential is to be set to zero for the given
     * atoms
     */
    public boolean isZero(CoordinatePair cPair) {
  		return isZero(cPair.r2());
    }
    
    public abstract double getRange();
        
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
    public abstract Potential0Lrc makeLrcPotential(Space space, Potential2 potential);
    

    ///************** end of methods for PotentialTruncation ***************
    
    public static final PotentialTruncation NULL = new Null();
    /**
     * No-op version of PotentialTruncation that performs no truncation at all.
     */
     private static final class Null extends PotentialTruncation {
        public boolean isZero(double r2) {return false;}
        public double uTransform(double r2, double untruncatedValue) {return untruncatedValue;}
        public double duTransform(double r2, double untruncatedValue) {return untruncatedValue;}
        public double d2uTransform(double r2, double untruncatedValue) {return untruncatedValue;}
        public Potential0Lrc makeLrcPotential(Space space, Potential2 potential) {return null;}
        public double getRange() {return Double.MAX_VALUE;}
     }//end of Null
}//end of PotentialTruncation