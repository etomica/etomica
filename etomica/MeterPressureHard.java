//this class includes a main method to demonstrate its use
package etomica;
import etomica.units.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.Frame;

/**
 * Meter for the pressure of a hard potential.
 * Performs sum of collision virial over all collisions, dividing by 
 * appropriate terms to obtain the pressure.
 */
public final class MeterPressureHard extends Meter implements IntegratorHard.CollisionListener, EtomicaElement {
        
    private double virialSum = 0.0;
    private double t0 = 0.0; //initialized in setPhaseIntegrator method
    private IntegratorHard integratorHard;
    private double value = 0.0; //holds meter value from one call to the next
    private final int D;
    
    public MeterPressureHard() {
        this(Simulation.instance);
    }
    public MeterPressureHard(Simulation sim) {
        super(sim);
        D = parentSimulation().space().D();
        setLabel("PV/Nk");
    }
        
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure measured via impulsive virial averaged over interatomic hard collisions");
        return info;
    }
        
    /**
     * Indicator that this meter returns quantity of dimension TEMPERATURE (note: returns PV/Nk, not P)
     */
    public Dimension getDimension() {return Dimension.TEMPERATURE;}

    /**
     * Returns P*V/N*kB = T - (virial sum)/(elapsed time)/(space dimension)/(number of atoms)
     * Virial sum and elapsed time apply to period since last call to this method.
     */
    public double currentValue() {
        double t = integratorHard.elapsedTime();
        if(t > t0) { // could have t==t0 if called twice before advancing integator
            double flux = -virialSum/((t-t0)*(double)(D*phase.atomCount()));   //divide by time interval
            value = integratorHard.temperature() + flux;
            t0 = t;
            virialSum = 0.0;
        }
        return value;
    }
    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    public void collisionAction(AtomPair pair, Potential.Hard p) {
        virialSum += p.lastCollisionVirial();
    }
    
    /**
     * Invokes superclass method and registers meter as a collisionListener to integrator.
     * Performs only superclass method if integrator is not an instance of IntegratorHard.
     */
	protected void setPhaseIntegrator(Integrator newIntegrator) {
	    super.setPhaseIntegrator(newIntegrator);
	    if(newIntegrator == null) return;
	    if(newIntegrator instanceof IntegratorHard) {
	        integratorHard = (IntegratorHard)newIntegrator;
	        integratorHard.addCollisionListener(this);
    	    t0 = integratorHard.elapsedTime();
	    }
	    else {  //should have an exception for this
	        System.out.println("Error in integrator type in MeterPressureHard");
	        System.exit(1);
	    }
	}
	
    /**
     * Method to demonstrate and test the use of this class.  
     * Pressure is measured in a hard-disk MD simulation.
     */
    public static void main(String[] args) {
        Frame f = new Frame();   //create a window
        f.setSize(600,350);
        Simulation.makeSimpleSimulation();  
        
        //here's the part unique to this class
        Phase phase = Simulation.instance.phase(0);
        Integrator integrator = Simulation.instance.integrator(0);
        integrator.setIsothermal(true);
        integrator.setTemperature(Kelvin.UNIT.toSim(300.));
        //make the meter and register it with the integrator
        MeterPressureHard pressureMeter = new MeterPressureHard();
        //Meter must be registered as collision listener and as interval listener to integrator
        //This is completed by the setPhase method
        //It is not be good to register the same listener multiple times, since addIntervalListener list does not prohibit redundant entries (addCollisionListener however does not muliply register the same listener)
           // done by setPhase:  ((IntegratorHard)integrator).addCollisionListener(pressureMeter);
           // done by setPhase:  integrator.addIntervalListener(pressureMeter);
           
        //set the phase where the meter performs its measurement and register as listener to phase's integrator
        //this call is commented out here since the setPhase call is performed by the default (but not by the Basic) elementCoordinator
        //note that there is no harm in calling setPhase multiple times with the same phase
           // pressureMeter.setPhase(phase);
           
        //display the meter
        DisplayBox box = new DisplayBox();
        box.setMeter(pressureMeter);
        //end of unique part
 
        Simulation.instance.elementCoordinator.go();
        phase.firstSpecies().setNMolecules(60);
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(WindowEvent e) {System.exit(0);}
        });
    }
    
}
