package etomica.data.meter;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHard;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.units.Pressure;

/**
 * Meter for the pressure (given as the compressibility factor) of a hard potential.
 * Performs sum of collision virial over all collisions, and manipulates value
 * to obtain the compressibility factor, PV/NkT.
 *
 * @author David Kofke
 */
public class MeterPressureHard extends DataSourceScalar implements
                                                IntegratorHard.CollisionListener,
                                                DataSourceCollisional {
    
    public MeterPressureHard(Space space) {
        super("Pressure", Pressure.dimension(space.D()));
    }
        
    /**
     * Returns P = (NT - (virial sum)/((elapsed time)*T*(space dimension)))/V
     * Virial sum and elapsed time apply to period since last call to this method.
     */
    public double getDataAsScalar() {
        if (integratorHard == null) throw new IllegalStateException("must call setIntegrator before using meter");
        if (!integratorHard.isIsothermal()) {
            throw new IllegalStateException("Integrator must be isothermal");
        }
        Phase phase = integratorHard.getPhase();
        double currentTime = integratorHard.getCurrentTime();
        double elapsedTime = currentTime - lastTime;
        if(elapsedTime == 0.0) return Double.NaN;
        double value = (integratorHard.getTemperature()*phase.atomCount() - virialSum/(phase.getSpace().D()*elapsedTime)) / 
                        phase.getBoundary().volume();

        virialSum = 0.0;
        lastTime = currentTime;
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
	public void setIntegrator(IntegratorHard newIntegrator) {
		if(newIntegrator == integratorHard) return;
		if(integratorHard != null) {
            integratorHard.removeCollisionListener(this);
        }
        integratorHard = newIntegrator;
	    if(newIntegrator != null) {
            integratorHard.addCollisionListener(this);
            lastTime = integratorHard.getCurrentTime();
        }
        virialSum = 0;
	}
    
    public IntegratorHard getIntegrator() {
        return integratorHard;
    }
    
    private static final long serialVersionUID = 1L;
    protected double virialSum;
    protected IntegratorHard integratorHard;
    protected double lastTime;
}
