package etomica;

/**
 * Harmonic Well interatomic potential.
 * Spherically symmetric potential of the form u(r) = 0.5*springConstant*(r)^2
 * where springConstant describes the strength of the pair interaction
 */
public class PotentialHarmonic extends Potential implements Potential.Soft {

    private Space.Vector force;
    private double w = 100.0;// Spring constant gives a measure of the strength of harmonic interaction
    
    public PotentialHarmonic(double w) {
        this(Simulation.instance, w);
    }
    public PotentialHarmonic(Simulation sim, double w) {
        super(sim);
        setW(w);
        force = sim.space().makeVector();
    }
   
   /**
    * Always returns false
    */
    public boolean overlap(AtomPair pair) {return false;}  //might want to change this

   /**
    * Energy of the given pair.
    */
    public double energy(AtomPair pair) {
        return 0.5*w*pair.r2();
    }
   /**
    * Force that atom2 exerts on atom1
    */
    public Space.Vector force(AtomPair pair) {
        force.E(pair.dr());
        force.TE(-w);
        return force;
    }
   
    /**
     * Virial, defined r*du/dr.
     * Used to compute the pressure
     */
    public double virial(AtomPair pair) {  //not carefully checked for correctness
        return -w*pair.r2();
    }
    
    /**
     * Returns the total (extensive) long-range correction to the energy
     * Input are the number of atoms of each type (i.e., x1*N and x2*N, where x1 and x2
     * are the mole fractions of each species interacting according to this potential), and the volume
     */
    public double energyLRC(int n1, int n2, double V) {
        return 0.0;
    }
    /**
     * Returns the total long-range correction to the pressure
     * Input are the number of atoms of each type (i.e., x1*N and x2*N, where x1 and x2
     * are the mole fractions of each species interacting according to this potential), and the volume
     */
    public double pressureLRC(int n1, int n2, double V) {
        return 0.0;
    }
    
    /**
     * Accessor method for harmonic energy parameter
     */
    public double getW() {return w;}
    /**
     * Accessor method for harmonic energy parameter
     */
    public void setW(double factor) {
        w = factor;
    }
    
}
  