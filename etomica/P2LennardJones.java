package etomica;
import etomica.units.Dimension;

/**
 * Lennard-Jones interatomic potential.
 * Spherically symmetric potential of the form u(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
 * where epsilon describes the strength of the pair interaction, and sigma is the atom size parameter
 */
public class P2LennardJones extends Potential2SoftSpherical implements EtomicaElement {

    public PotentialTruncation truncation;
    
  public String getVersion() {return "PotentialLJ:01.07.03/"+Potential.VERSION;}

    public P2LennardJones() {
        this(Simulation.instance, Default.ATOM_SIZE, Default.POTENTIAL_WELL, Default.POTENTIAL_CUTOFF);
    }
    public P2LennardJones(double sigma, double epsilon, double cutoff) {
        this(Simulation.instance, sigma, epsilon, cutoff);
    }
    public P2LennardJones(Simulation sim, double sigma, double epsilon, double cutoff) {
        super(sim);
        setSigma(sigma);
        setEpsilon(epsilon);
        setCutoff(cutoff);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple Lennard-Jones potential");
        return info;
    }

    public double u(double r2) {
        if(truncation(r2).isZero) {return 0.0;}
        if(r2 > cutoffRadiusSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return epsilon4*s6*(s6 - 1.0);
        }
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {  //not carefully checked for correctness
        if(r2 > cutoffRadiusSquared) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return -epsilon48*s6*(s6 - 0.5);
        }
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
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
    public double uLRC() {return uLRC;}

 /*   public double uLRC(int n1, int n2, double V) {
        return n1*n2*uLRC/V;
    }*/
    /**
     * Returns the total long-range correction to the pressure
     * Input are the number of atoms of each type (i.e., x1*N and x2*N, where x1 and x2
     * are the mole fractions of each species interacting according to this potential), and the volume
     */
    public double duLRC() {return duLRC;}
/*    public double pressureLRC(int n1, int n2, double V) {
        return n1*n2*pLRC/(V*V);
    }*/

    //not implemented
    public double d2uLRC() {return 0.0;}
    
    
    
    /**
     * Accessor method for Lennard-Jones size parameter.
     */
    public double getSigma() {return sigma;}
    /**
     * Mutator method for Lennard-Jones size parameter.
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
        setCutoff(cutoff);
    }
    public Dimension getSigmaDimension() {return Dimension.LENGTH;}

    /**
     * Accessor method for potential cutoff class.
     */
    public PotentialTruncation getTruncation() {return truncation;}
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
        uLRC = 2.0*epsilon*sigmaD*A*(rc12/(12.-D) - rc6/(6.-D))/rcD;  //complete LRC is obtained by multiplying by N1*N2/rho
        duLRC = 2.0*epsilon*sigmaD*A*(12.*rc12/(12.-D) - 6.*rc6/(6.-D))/(D*rcD);  //complete LRC is obtained by multiplying by N1*N2/rho
        d2uLRC = 0.0;  //complete LRC is obtained by multiplying by N1*N2/rho
    }
    
    private double sigma, sigmaSquared;
    private double cutoffRadius, cutoffRadiusSquared;
    private double epsilon, cutoff;
    private double epsilon4, epsilon48, epsilon144;
    private double uLRC, duLRC, d2uLRC;  //multipliers for long-range corrections
}
  