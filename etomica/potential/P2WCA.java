package etomica.potential;
import etomica.Default;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Space;
import etomica.potential.Potential2SoftSpherical;
import etomica.units.Dimension;

/**
 * Weeks-Chandler-Andersen potential.  Obtained by truncating the Lennard-Jones
 * potential at the separation where it has its minimum, and shifting it upwards so the
 * minimum is at zero energy.  The resulting potential has a value and first
 * derivative that is continuous at the point of truncation.
 *
 * @author David Kofke (edited by Eric C. Cichowski and Todd Schmidt)
 */
public class P2WCA extends Potential2SoftSpherical implements EtomicaElement {

    /**
     * Constructs potential using default sigma and epsilon given by Default class.
     */
    public P2WCA(Space space) {
        this(space, Default.ATOM_SIZE, Default.POTENTIAL_WELL);
    }
    
    public P2WCA(Space space, double sigma, double epsilon) {
        super(space);
        setSigma(sigma);
        setEpsilon(epsilon);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Weeks-Chandler-Andersen potential");
        return info;
    }
    
    /**
     * Returns the range of the potential, which is the point of truncation.  This
     * is equal to 2^(1/6) * sigma.
     */
    public double getRange() {
        return range;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        if(r2 < rangeSquared) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            return epsilon4*s6*(s6 - 1.0) + epsilon;
        }
        else return 0.0;
     }


    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        if(r2 < rangeSquared) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            return -epsilon48*s6*(s6 - 0.5);
        }
        else return 0.0;
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        if(r2 < rangeSquared) {
            double s2 = sigmaSquared/r2;
            s6 = s2*s2*s2;
            return epsilon624*s6*(s6 - _168div624);
        }
        else return 0.0;
      }
            
    /**
     *  Returns zero, since there is no long-range correction.
     */
    public double uInt(double rC) {
        return 0.0;
    }

    /**
     * Accessor method for the size parameter.
     */
    public double getSigma() {return sigma;}
    /**
     * Mutator method for Lennard-Jones size parameter.
     * Does not adjust potential cutoff for change in sigma.
     */
    public final void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
		range = sigma*Math.pow(2,1./6.);
		rangeSquared = range*range;

    }
    public Dimension getSigmaDimension() {return Dimension.LENGTH;}
    
    /**
     * Accessor method for the energy parameter
     */
    public double getEpsilon() {return epsilon;}
    
    /**
     * Mutator method for the energy parameter
     */
    public final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = eps*4.0;
        epsilon48 = eps*48.0;
        epsilon624 = eps*624.0;
    }
    public Dimension getEpsilonDimension() {return Dimension.ENERGY;}
   
    private double sigma, sigmaSquared, range, rangeSquared;
    private double epsilon;
    private double epsilon4, epsilon48, epsilon624;
    private static final double _168div624 = 168./624.;
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
	    Potential2 potential11 = new P2LJWCA(1.0, 100.);//sigma, epsilon
	    Potential2 potential12 = new P2LJWCA(2.0, 600.);
	    Potential2 potential22 = new P2LJWCA(3.0, 300.);
	    
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
    
}//end of P2LJWCA
  