package etomica.potential;

import etomica.EtomicaElement;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Finite elastic nonlinear extensible (FENE) spring potential.
 * Intramolecular potential for modeling polymer chains.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/31/02 (DAK) new
  */

public class P2Fene extends Potential2SoftSpherical implements EtomicaElement {

    private double r0, r02, h, prefactor;
    
    public P2Fene(Simulation sim) {
        this(sim.space, 1.50*sim.getDefaults().atomSize, 
                30.0*sim.getDefaults().potentialWell/Math.pow(sim.getDefaults().atomSize,2)/*1.50*Default.atomSize, 250.0*/);
    }
    public P2Fene(Space space, double r0, double amplitude) {
        super(space);
        setMaximumSeparation(r0);
        setAmplitude(amplitude);
    }

    public double u(double r2) {
        return (r2 < r02) ? prefactor * Math.log(1 - r2/r02) : Double.POSITIVE_INFINITY;
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        return h*r2*r02/(r02 - r2);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        double d = (r02 - r2);
        return h*r2*r02*(r02 + r2)/(d*d);
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        return 0.0;
    }
    
    public double getAmplitude() {return h;}
    public void setAmplitude(double H) {
        this.h = H;
        prefactor = -0.5*h*r02-Math.log(r02);
    }
    
    /**
     * Not implemented correctly.  
     * Should be energy/length^3.
     */
    public Dimension getAmplitudeDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION,Length.DIMENSION},new double[]{1,-3});
    }

    /**
     * Accessor method for harmonic energy parameter
     */
    public double getMaximumSeparation() {return r0;}
    /**
     * Accessor method for harmonic energy parameter
     */
    public void setMaximumSeparation(double r0) {
        this.r0 = r0;
        r02 = r0*r0;
        prefactor = -0.5*h*r02;
    }
    
    public Dimension getMaximumSeparationDimension() {
        return Length.DIMENSION;
    }
    
}
