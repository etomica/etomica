package etomica.data.meter;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.data.DataSourceCountTime;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHard;
import etomica.space.Space;
import etomica.units.Dimension;

/**
 * Meter for the pressure (given as the compressibility factor) of a hard potential.
 * Performs sum of collision virial over all collisions, and manipulates value
 * to obtain the compressibility factor, PV/NkT.
 *
 * @author David Kofke
 */
public class MeterPressureHard extends DataSourceScalar implements
                                                IntegratorHard.CollisionListener,
                                                MeterCollisional,
                                                EtomicaElement {
    
    public MeterPressureHard(Space space, IntegratorHard integrator) {
        super("Pressure", Dimension.pressure(space.D()));
        timer = new DataSourceCountTime();
        setIntegrator(integrator);
    }
        
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Compressibility factor measured via impulsive virial averaged over interatomic hard collisions");
        return info;
    }
        
    /**
     * Returns P = (NT - (virial sum)/((elapsed time)*T*(space dimension)))/V
     * Virial sum and elapsed time apply to period since last call to this method.
     */
    //TODO consider how to ensure timer is advanced before this method is invoked
    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        if (integratorHard.getPhase()[0] != phase) {
            throw new IllegalArgumentException("Integrator must act on meter's phase");
        }
        else if (!integratorHard.isIsothermal()) {
            throw new IllegalStateException("Integrator must be isothermal");
        }
        double elapsedTime = timer.getDataAsScalar();
        if(elapsedTime == 0.0) return Double.NaN;
        double value = (integratorHard.getTemperature()*phase.atomCount() - virialSum/(phase.space().D()*elapsedTime)) / phase.boundary().volume();

        virialSum = 0.0;
        timer.reset();
        return value;
    }
    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        virialSum += agent.collisionPotential.lastCollisionVirial();
    }
    
    /**
     * Implementation of Meter.MeterCollisional interface.  Returns -(collision virial).
     * Suitable for tabulation of PV
     */
	public double collisionValue(IntegratorHard.Agent agent) {
	    return -agent.collisionPotential.lastCollisionVirial();
	}

    /**
     * Registers meter as a collisionListener to the integrator, and sets up
     * a DataSourceTimer to keep track of elapsed time of integrator.
     */
	protected void setIntegrator(IntegratorHard newIntegrator) {
		if(newIntegrator == integratorHard) return;
		if(integratorHard != null) {
            integratorHard.removeCollisionListener(this);
            integratorHard.removeListener(timer);
        }
        integratorHard = newIntegrator;
	    if(newIntegrator != null) {
            newIntegrator.addListener(timer);
            integratorHard.addCollisionListener(this);
        }
	}
    
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    protected double virialSum = 0.0;
    private IntegratorHard integratorHard;
    protected final DataSourceCountTime timer;

	
    /**
     * Method to demonstrate and test the use of this class.  
     * Pressure is measured in a hard-sphere MD simulation.
     */
 /*   public static void main(String[] args) {
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        
        //here's the part unique to this class
        sim.integrator.setIsothermal(true);
        sim.integrator.setTemperature(Kelvin.UNIT.toSim(300.));
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
           
        MeterProfileHard profile = new MeterProfileHard();
        profile.setActive(true);
        profile.setX(0,30,50);
        profile.setMeter(pressureMeter);
        profile.setMeter(((MeterScalar.MeterScalarAtomic)new MeterTemperature()));
        etomica.graphics.DisplayPlot plot = new etomica.graphics.DisplayPlot();
           
        //display the meter
        etomica.graphics.DisplayBox box = new etomica.graphics.DisplayBox();
        box.setMeter(pressureMeter);
        
        etomica.action.MeterReset resetAction = new etomica.action.MeterReset();
        resetAction.setMeters(new MeterAbstract[] {pressureMeter, profile});
        etomica.graphics.DeviceButton resetButton = new etomica.graphics.DeviceButton(resetAction);
        //end of unique part
 
        Simulation.instance.elementCoordinator.go();
        plot.setDataSources(profile);
        sim.phase.firstSpecies().setNMolecules(60);
        
        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
    }//end of main
    */
}
