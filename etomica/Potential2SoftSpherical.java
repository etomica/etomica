package etomica;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 */
public interface Potential2SoftSpherical {
        
   /**
    * The pair energy u(r^2).
    * @param the square of the distance between the particles.
    */
    public double energy(double r2);
        
   /**
    * The derivative of the pair energy, times the separation r: r du/dr.
    */
    public double rDuDr(double r2);
        
   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double r2D2uDr2(double r2);
        
   /**
    * Long-range contribution to the energy, calculated assuming a uniform distribution
    * of other particles from the cutoff to infinity.
    */
    public double uLRC();
        
   /**
    * Long-range contribution to the energy first derivative.
    */
    public double rDuDrLRC();
        
   /**
    * Long-range contribution to the energy second derivative.
    */
    public double r2D2uDr2LRC();
    
}//end of Potential2SoftSpherical
