package etomica;

/**
 * Class that defines whether and how the interatomic potential is truncated.
 *
 * @author David Kofke
 */
 
public abstract class PotentialTruncation {
    
    /**
     * Returns true is the truncation makes the potential zero.
     */
    public abstract boolean isZero(double r2);
    /**
     * Transforms the given value of the energy for the truncation applied at
     * the given (square) distance.
     */
    public abstract double uTransform(double r2, double untruncatedValue);

    ///************** end of methods fo PotentialTruncation ***************
    
    public static final PotentialTruncation NULL = new Null();
    /**
     * No-op version of PotentialTruncation that performs no truncation at all.
     */
     private static final class Null extends PotentialTruncation {
        public boolean isZero(double r2) {return false;}
        public double uTransform(double r2, double untruncatedValue) {
            return untruncatedValue;
        }
     }//end of Null
        
    /**
     * Simple truncation of the potential as an adjustable cutoff separation.
     * The energy is unaffected for separations less than the truncation distance,
     * and it is set to zero beyond this distance.
     */
    public static final class Simple extends PotentialTruncation {
        
        double rCutoff, r2Cutoff;
        
        public Simple() {this(Default.POTENTIAL_CUTOFF);}
        public Simple(double rCutoff) {setCutoff(rCutoff);}
        public boolean isZero(double r2) {return r2 > r2Cutoff;}
        public double uTransform(double r2, double untruncatedValue) {
            return (r2 > r2Cutoff) ? 0.0 : untruncatedValue;
        }
        public final void setCutoff(double rCut) {
            rCutoff = rCut;
            r2Cutoff = rCut*rCut;
        }
        public double getCutoff() {return rCutoff;}
        
    }//end of Simple
}//end of PotentialTruncation