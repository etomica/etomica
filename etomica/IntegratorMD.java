package etomica;

import etomica.units.*;
/**
 * Superclass of all molecular-dynamics integrators.
 * Extends the Integrator class by adding methods and inner classes that 
 * set and manipulate the time step and elapsed simulation time.
 */
public abstract class IntegratorMD extends Integrator {
    
    /**
     * Elementary time step for the MD simulation
     */
    protected double timeStep;
    /**
     * Amount of simulation time elapsed since the start of the integration
     */
    protected double elapsedTime;
    
    /**
     * Zero-point counter used to manage the elapsedTime method
     */
    private int stepCount0 = 0;
    /**
     * Zero-point counter used to manage the elapsedTime method
     */
    private double t0 = 0.0;
    
    public IntegratorMD(PotentialMaster potentialMaster) {
        super(potentialMaster);
        setTimeStep(Default.TIME_STEP);
    }
    
    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return false;
//        potential.set(firstPhase);  //assumes integrator has exclusive access to the potential hierarchy
        return true;
    }

    /**
     * The simulation time elapsed since the start of the integration.
     * Cannot be reset to zero.
     */
//    public double elapsedTime() {return t0 + (stepCount-stepCount0)*timeStep;}
    
    /**
     * Zeros step counter and reference time (such that elapsed time becomes zero).
     * Often overridden in subclasses, so that elapsed time does not get zeroed
     * with iterator reset.
     */
    protected void doReset() {
//        stepCount = 0;
//        stepCount0 = 0;
        t0 = 0;
    }
        
    /**
     * Sets integration time step.
     * Updates zero-point counters used internally to manage the elapsedTime method
     */
    public void setTimeStep(double t) {
        timeStep = t;
//        t0 = elapsedTime();
//        stepCount0 = stepCount;
    }
    public final double getTimeStep() {return timeStep;}
    public final double timeStep() {return timeStep;}
    public Dimension getTimeStepDimension() {return Dimension.TIME;}
    
}//end of IntegratorMD
    