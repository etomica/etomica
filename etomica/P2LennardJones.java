package etomica;
import etomica.units.Dimension;

/**
 * Lennard-Jones interatomic potential.
 * Spherically symmetric potential of the form u(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
 * where epsilon describes the strength of the pair interaction, 
 * and sigma is the atom size parameter.
 *
 * @author David Kofke
 */
public final class P2LennardJones extends Potential2SoftSpherical implements EtomicaElement {

  public final PotentialTruncation truncation;
    
  public String getVersion() {return "PotentialLJ:01.07.05/"+Potential2SoftSpherical.VERSION;}

    public P2LennardJones() {
        this(Default.ATOM_SIZE, Default.POTENTIAL_WELL);
    }
    public P2LennardJones(double sigma, double epsilon) {
        super(Simulation.instance);
        setSigma(sigma);
        setEpsilon(epsilon);
        truncation = new PotentialTruncationSimple(this, Default.POTENTIAL_CUTOFF_FACTOR * sigma);
    }
    public P2LennardJones(Simulation sim, double sigma, double epsilon,
                            PotentialTruncation trunc) {
        super(sim);
        setSigma(sigma);
        setEpsilon(epsilon);
        truncation = trunc;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple Lennard-Jones potential");
        return info;
    }

    public double u(double r2) {
        if(truncation.isZero(r2)) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            //may need to restructure to remove overhead of method call
            return truncation.uTransform(r2, epsilon4*s6*(s6 - 1.0));
        }
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        if(truncation.isZero(r2)) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return truncation.duTransform(r2, -epsilon48*s6*(s6 - 0.5));
        }
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        if(truncation.isZero(r2)) {return 0.0;}
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            return truncation.d2uTransform(r2, epsilon624*s6*(s6 - _168div624));
        }
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        if(rC != rCLast) { //recompute if something changed, otherwise used saved value
            rCLast = rC;
            double A = parentSimulation().space().sphereArea(1.0);  //multiplier for differential surface element
            int D = parentSimulation().space().D();                 //spatial dimension
            double sigmaD = 1.0;  //will be sigma^D
            double rcD = 1.0;     //will be (sigam/rc)^D
            double rc = sigma/rC;
            for(int i=D; i>0; i--) {
                sigmaD *= sigma;
                rcD *= rc;
            }
            double rc3 = rc*rc*rc;
            double rc6 = rc3*rc3;
            double rc12 = rc6*rc6;
            uInt = 4.0*epsilon*sigmaD*A*(rc12/(12.-D) - rc6/(6.-D))/rcD;  //complete LRC is obtained by multiplying by N1*N2/rho
        }
        return uInt;
    }

    /**
     * Accessor method for potential cutoff class.
     */
    public PotentialTruncation getTruncation() {return truncation;}
    
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
        rCLast = 0.0;
    }
    public Dimension getSigmaDimension() {return Dimension.LENGTH;}
    
    /**
     * Accessor method for Lennard-Jones energy parameter
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Mutator method for Lennard-Jones energy parameter
     */
    public final void setEpsilon(double eps) {
        uInt *= eps/epsilon;
        epsilon = eps;
        epsilon4 = eps*4.0;
        epsilon48 = eps*48.0;
        epsilon624 = eps*624.0;
    }
    public Dimension getEpsilonDimension() {return Dimension.ENERGY;}
   
    private double sigma, sigmaSquared;
    private double epsilon;
    private double epsilon4, epsilon48, epsilon624;
    private static final double _168div624 = 168./624.;
    private double uInt, rCLast;  
}//end of P2LennardJones
  