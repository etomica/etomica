package etomica;
import etomica.units.Dimension;

/**
 * Lennard-Jones interatomic potential.
 * Spherically symmetric potential of the form u(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
 * where epsilon describes the strength of the pair interaction, and sigma is the atom size parameter
 */
public class PotentialLJ extends Potential2 implements Potential.Soft, EtomicaElement {

  public String getVersion() {return "PotentialLJ:01.05.25/"+Potential.VERSION;}

    private double sigma, sigmaSquared;
    private double cutoffRadius, cutoffRadiusSquared;
    private double epsilon, cutoff;
    private double epsilon4, epsilon48, epsilon144;
    private double eLRC, pLRC;  //multipliers for long-range correction to energy and pressure, resp.
    private final Space.Vector force;

    public PotentialLJ() {
        this(Simulation.instance, Default.ATOM_SIZE, Default.POTENTIAL_WELL, Default.POTENTIAL_CUTOFF);
    }
    public PotentialLJ(double sigma, double epsilon, double cutoff) {
        this(Simulation.instance, sigma, epsilon, cutoff);
    }
    public PotentialLJ(Simulation sim, double sigma, double epsilon, double cutoff) {
        super(sim);
        setSigma(sigma);
        setCutoff(cutoff);
        setEpsilon(epsilon);
        force = sim.space().makeVector();
        calculateLRC();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple Lennard-Jones potential");
        return info;
    }

   /**
    * Always returns false
    */
    public boolean overlap(AtomPair pair) {return false;}  //might want to change this

   /**
    * Energy of the given pair.
    * Returns zero if greater than the cutoff separation
    */
    public double energy(AtomPair pair) {
        double r2 = pair.r2();
        if(r2 > cutoffRadiusSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return epsilon4*s6*(s6 - 1.0);
        }
    }
   /**
    * Force that atom2 exerts on atom1
    */
    public Space.Vector force(AtomPair pair) {
        double r2 = pair.r2();
        if(r2 > cutoffRadiusSquared) {force.E(0.0);}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            force.E(pair.dr());
            force.TE(-epsilon48*s6*(s6-0.5)/r2);
        }
        return force;
    }            
   
    /**
     * Virial, defined r*du/dr.
     * Used to compute the pressure
     */
    public double virial(AtomPair pair) {  //not carefully checked for correctness
        double r2 = pair.r2();
        if(r2 > cutoffRadiusSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return -epsilon48*s6*(s6 - 0.5);
        }
    }
    
    /**
     * Hypervirial, defined r*du/dr + r^2*d2u/dr2.
     * Used to compute the pressure
     */
    public double hyperVirial(AtomPair pair) {  //not carefully checked for correctness
        double r2 = pair.r2();
        if(r2 > cutoffRadiusSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return epsilon144*s6*(4.*s6 - 1.0);
        }
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
     * Accessor method for Lennard-Jones size parameter
     */
    public double getSigma() {return sigma;}
    /**
     * Accessor method for Lennard-Jones size parameter
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
        setCutoff(cutoff);
    }
    public Dimension getSigmaDimension() {return Dimension.LENGTH;}

    /**
     * Accessor method for Lennard-Jones cutoff distance; divided by sigma
     * @return cutoff distance, divided by size parameter (sigma)
     */
    public double getCutoff() {return cutoff;}
    /**
     * Accessor method for Lennard-Jones cutoff distance; divided by sigma
     * @param rc cutoff distance, divided by size parameter (sigma)
     */
    public final void setCutoff(double rc) {  
        cutoff = rc;
        cutoffRadius = sigma*cutoff;
        cutoffRadiusSquared = cutoffRadius*cutoffRadius;
        calculateLRC();
    }
    public Dimension getCutoffDimension() {return Dimension.LENGTH;}
    
    /**
     * Accessor method for Lennard-Jones energy parameter
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Accessor method for Lennard-Jones energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = eps*4.0;
        epsilon48 = eps*48.0;
        epsilon144 = eps*144.0;
        calculateLRC();
    }
    public Dimension getEpsilonDimension() {return Dimension.ENERGY;}
   
    //Calculates basic parameters for returning long-range correction to energy and virial
    private void calculateLRC() {
        double A = parentSimulation().space().sphereArea(1.0);  //multiplier for differential surface element
        int D = parentSimulation().space().D();                         //spatial dimension
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
    }
    
}
  