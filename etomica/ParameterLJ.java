package etomica;
import etomica.units.Dimension;

/**
 * Parameter having two fields, one of dimension length, and one of dimension
 * energy.  These are the fields needed to define a Lennard-Jones potential.
 *
 * @author David Kofke
 */
 
 public class ParameterLJ extends Parameter implements Parameter.Size, Parameter.Energy {
    
    /**
     * Accessor method for the size parameter.
     */
    public double getSigma() {return sigma;}
    
    /**
     * Mutator method for size parameter.
     */
    public final void setSigma(double s) {sigma = s;}
    
    /**
     * Reports sigma as having dimensions of length.
     */
    public Dimension getSigmaDimension() {return Dimension.LENGTH;}
    
    /**
     * Accessor method for energy parameter
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Mutator method for energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
    }
    /**
     * Reports epsilon as having dimensions of energy.
     */
    public Dimension getEpsilonDimension() {return Dimension.ENERGY;}
   
    private double sigma, epsilon;
 }   