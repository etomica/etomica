package etomica;

import etomica.units.*;
/**
 * Superclass of all molecular-dynamics integrators.
 * Extends the Integrator class by adding methods and inner classes that 
 * set and manipulate the time step and elapsed simulation time.
 */
public abstract class IntegratorMD extends Integrator {
    
    public static String VERSION = "IntegratorMD:01.03.01.0/"+Integrator.VERSION;
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
    
    public IntegratorMD(SimulationElement parent) {
        super(parent);
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
    public double elapsedTime() {return t0 + (stepCount-stepCount0)*timeStep;}
    
    /**
     * Zeros step counter and reference time (such that elapsed time becomes zero).
     * Often overridden in subclasses, so that elapsed time does not get zeroed
     * with iterator reset.
     */
    protected void doReset() {
        stepCount = 0;
        stepCount0 = 0;
        t0 = 0;
    }
        
    /**
     * Sets integration time step.
     * Updates zero-point counters used internally to manage the elapsedTime method
     */
    public void setTimeStep(double t) {
        timeStep = t;
        t0 = elapsedTime();
        stepCount0 = stepCount;
    }
    public final double getTimeStep() {return timeStep;}
    public final double timeStep() {return timeStep;}
    public Dimension getTimeStepDimension() {return Dimension.TIME;}
    
    /**
     * Returns a Meter object that gives the elapsed simulation time.
     * @deprecated use MeterTime
     */
    public ChronoMeter chronoMeter() {
        return new ChronoMeter();
    }
    
    /**
     * Class for measuring elapsed simulation time
     *
     * @deprecated  Use MeterTime instead.
     */
    public class ChronoMeter extends MeterScalar {
        
        double time0;
        
        /**
         * Constructor for a time-meter associated with this IntegratorMD class.
         */
        ChronoMeter() {
            super(IntegratorMD.this.simulation());
            setActive(false);
            this.reset();
            setLabel("Elapsed time");
        }
        
        /**
         * Returns the integrator for which this meter keeps time.
         */
        public Integrator integrator() {return IntegratorMD.this;}
        
        /**
         * Indicates units (time) of measured property.
         */
        public Dimension getDimension() {return Dimension.TIME;}
        
        /**
         * Sets the meter's current time as the new zero for subsequent values.
         */
        public void reset() {time0 = elapsedTime();}
        
        /**
         * Overrides parent class method to cause no action to be performed in response 
         * to interval event.
         */
        public void intervalAction(Integrator.IntervalEvent evt) {}
        /**
         * Returns the simulation time elapsed since the instantiation of
         * this meter, or since the last call to reset().
         */
        public double currentValue() {return elapsedTime() - time0;}
        /**
         * Same as currentValue().
         */
    	public double mostRecent() {return currentValue();}
    }//end of ChronoMeter
}//end of IntegratorMD
    