package etomica;

import etomica.electrostatics.*;

/**
 * Potential for the primitive model of an electrolyte.
 * The primitive model is a simple Coulombic charge with a hard-sphere repulsion.
 * Charges are ascribed to the Atoms via their AtomType.electroType field
 * Hard-sphere diameter is defined by this potential
 */
public class PotentialPrimitive extends Potential implements Potential.Soft {

    private double sigma, sigmaSquared;
    private double cutoffRadius, cutoffRadiusSquared;
    private double cutoff;
    private double eLRC, pLRC;  //multipliers for long-range correction to energy and pressure, resp.
    private final Space.Vector force;

    public PotentialPrimitive() {
        this(Simulation.instance, Default.ATOM_SIZE, Default.POTENTIAL_CUTOFF);
    }
    public PotentialPrimitive(double sigma, double cutoff) {
        this(Simulation.instance, sigma, cutoff);
    }
    public PotentialPrimitive(Simulation sim, double sigma, double cutoff) {
        super(sim);
        setSigma(sigma);
        setCutoff(cutoff);
        force = sim.space().makeVector();
    }
 
   /**
    * Returns true if separation is less than sigma, false otherwise.
    */
    public boolean overlap(AtomPair pair) {return pair.r2() < sigmaSquared;}
   /**
    * Returns primitive-model energy.
    * Return infinity if overlap is true, and zero if separation is greater than cutoff.
    */
    public double energy(AtomPair pair) {
        double r2 = pair.r2();
        if(r2 > cutoffRadiusSquared) {return 0.0;}
        else if(r2 < sigmaSquared) {return Double.MAX_VALUE;}
        else {
            double z1 = ((Monopole)pair.atom1.type.electroType()).getZ();
            double z2 = ((Monopole)pair.atom2.type.electroType()).getZ();
//            return z1*z2/Math.sqrt(r2);
            return -0.5*z1*z2*Math.log(r2);
        }
    }
    /** 
     * Force that atom2 exerts on atom1
     * This method has not been checked and may be incorrect
     */
    public Space.Vector force(AtomPair pair) {
        double r2 = pair.r2();
        if(r2 > cutoffRadiusSquared) {force.E(0.0);}
        else if(r2 < sigmaSquared) {force.E(Double.MAX_VALUE);}
        else {
            double z1 = ((Monopole)pair.atom1.type.electroType()).getZ();
            double z2 = ((Monopole)pair.atom2.type.electroType()).getZ();
            double c = z1*z2/Math.sqrt(r2);
            force.E(pair.dr());
            force.TE(-c/r2);
        }
        return force;
    }            
     
     /**
      * Always returns zero. Not yet implemented.
      */
    public double virial(AtomPair pair) {  //not carefully checked for correctness
        return 0.0;
//        double r2 = pair.r2();
//        if(r2 > cutoffRadiusSquared) {return 0.0;}
//        else {
//            double s2 = sigmaSquared/r2;
//            double s6 = s2*s2*s2;
//            return -epsilon48*s6*(s6 - 0.5);
//        }
    }
    
    /**
     * Returns the total (extensive) long-range correction to the energy
     * Input are the number of atoms of each type (i.e., x1*N and x2*N, where x1 and x2
     * are the mole fractions of each species interacting according to this potential), and the volume
     */
    public double energyLRC(int n1, int n2, double V) {
        return n1*n2*eLRC/V;
    }
    /**
     * Returns the total long-range correction to the pressure
     * Input are the number of atoms of each type (i.e., x1*N and x2*N, where x1 and x2
     * are the mole fractions of each species interacting according to this potential), and the volume
     */
    public double pressureLRC(int n1, int n2, double V) {
        return n1*n2*pLRC/(V*V);
    }
    
    
    /**
     * Accessor method for the size of the repulsive core of the primitive model
     */
    public double getSigma() {return sigma;}
    /**
     * Accessor method for the size of the repulsive core of the primitive model
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
        setCutoff(cutoff);
    }

    /**
     * Accessor method for cutoff distance; divided by sigma
     * @return cutoff distance, divided by size parameter (sigma)
     */
    public double getCutoff() {return cutoff;}
    /**
     * Accessor method for cutoff distance; divided by sigma
     * @param rc cutoff distance, divided by size parameter (sigma)
     */
    public final void setCutoff(double rc) {  //argument is the cutoff radius divided by sigma
        cutoff = rc;
        cutoffRadius = sigma*cutoff;
        cutoffRadiusSquared = cutoffRadius*cutoffRadius;
        calculateLRC();
    }
    
    //Calculates basic parameters for returning long-range correction to energy and virial
    private void calculateLRC() {
/* from LJ       if(parentSimulation == null) return;
        double A = parentSimulation.space.sphereArea(1.0);  //multiplier for differential surface element
        int D = parentSimulation.D;                         //spatial dimension
        double sigmaD = 1.0;  //will be sigma^D
        double rcD = 1.0;     //will be (sigam/rc)^D
        double rc = sigma/cutoffRadius;
        for(int i=D; i>0; i--) {
            sigmaD *= sigma;
            rcD *= rc;
        }
        double rc3 = rc*rc*rc;
        double rc6 = rc3*rc3;
        double rc12 = rc6*rc6;
        eLRC = 2.0*epsilon*sigmaD*A*(rc12/(12.-D) - rc6/(6.-D))/rcD;  //complete LRC is obtained by multiplying by N1*N2/rho
        pLRC = 2.0*epsilon*sigmaD*A*(12.*rc12/(12.-D) - 6.*rc6/(6.-D))/(D*rcD);  //complete LRC is obtained by multiplying by N1*N2/rho
*/
    }
    
}
  