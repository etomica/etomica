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

    public P2LennardJones() {
        this(Default.ATOM_SIZE, Default.POTENTIAL_WELL);
    }
    public P2LennardJones(Space space) {
        super(space);
        setSigma(Default.ATOM_SIZE);
        setEpsilon(Default.POTENTIAL_WELL);
    }
    public P2LennardJones(double sigma, double epsilon) {
        super(Simulation.getDefault().space);//can't "this" with other constructor
        setSigma(sigma);
        setEpsilon(epsilon);
    }
    public P2LennardJones(Space space, double sigma, double epsilon,
                            PotentialTruncation trunc) {
        super(space, trunc);
        setSigma(sigma);
        setEpsilon(epsilon);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple Lennard-Jones potential");
        return info;
    }

    /**
     * The energy u.  No truncation is applied here; 
     * instead it is applied in the energy(AtomPair) method of Potential2SoftSpherical.
     */
    public double u(double r2) {
        if(r2 != r2Last) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            r2Last = r2;
        }
        return epsilon4*s6*(s6 - 1.0);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        if(r2 != r2Last) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            r2Last = r2;
        }
        return -epsilon48*s6*(s6 - 0.5);
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        if(r2 != r2Last) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            r2Last = r2;
        }
        return epsilon624*s6*(s6 - _168div624);
    }
            
    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) {
        if(rC != rCLast) { //recompute if something changed, otherwise used saved value
            rCLast = rC;
            double A = space.sphereArea(1.0);  //multiplier for differential surface element
            int D = space.D();                 //spatial dimension
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
     * Accessor method for Lennard-Jones size parameter.
     */
    public double getSigma() {return sigma;}
    /**
     * Mutator method for Lennard-Jones size parameter.
     * Does not adjust potential cutoff for change in sigma.
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
    private double r2Last = -1.0;
    private double s6;
    
    /**
     * main method to test and demonstrate use of potential.
     * Simulates a two-component mixture.
     */
/*    public static void main(String[] args) {
        
        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new Space2D());
        Simulation.instance = sim;
        
	    Phase phase = new Phase();
	    
	    //set up species
	    SpeciesSpheresMono species1 = new SpeciesSpheresMono();
	    SpeciesSpheresMono species2 = new SpeciesSpheresMono();
	    species1.setDiameter(1.0);
	    species2.setDiameter(3.0);
	    
	    //define potentials
	    Potential2 potential11 = new P2LennardJones(1.0, 100.);//sigma, epsilon
	    Potential2 potential12 = new P2LennardJones(2.0, 600.);
	    Potential2 potential22 = new P2LennardJones(3.0, 300.);
	    
	    //connect species and potentials
	    potential11.setSpecies(species1, species1);
	    potential12.setSpecies(species1, species2);
	    potential22.setSpecies(species2, species2);
	    
	    //control
	    Controller controller = new Controller();
	    IntegratorMD integrator = new IntegratorVelocityVerlet();
	    
	    //display
	    etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase();
	    etomica.graphics.DisplayTimer timer = new etomica.graphics.DisplayTimer(integrator);
	    timer.setUpdateInterval(10);
	    
	    display.setColorScheme(new etomica.graphics.ColorSchemeByType());
	    etomica.graphics.ColorSchemeByType.setColor(species1, java.awt.Color.blue);
		
        etomica.graphics.DeviceTrioControllerButton button = new etomica.graphics.DeviceTrioControllerButton();

        //meters
		MeterEnergy energy = new MeterEnergy();
		energy.setHistorying(true);
		energy.setActive(true);
		
		energy.getHistory().setNValues(500);
		
		etomica.graphics.DisplayPlot plot = new etomica.graphics.DisplayPlot();
		plot.setLabel("Energy");
		plot.setDataSource(energy.getHistory());
		
		integrator.setSleepPeriod(2);
		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
    } // */
    
}//end of P2LennardJones
  