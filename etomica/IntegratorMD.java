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
    
    public IntegratorMD(Simulation sim) {
        super(sim);
        setTimeStep(Default.TIME_STEP);
    }
    
    /**
     * The simulation time elapsed since the start of the integration.
     * Cannot be reset to zero.
     */
    public double elapsedTime() {return t0 + (stepCount-stepCount0)*timeStep;}
        
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
     */
    public ChronoMeter chronoMeter() {
        return new ChronoMeter();
    }
    
    /**
     * Class for measuring elapsed simulation time
     */
    public class ChronoMeter extends Meter {
        
        double t0;
        
        ChronoMeter() {
            super(IntegratorMD.this.parentSimulation());
            setActive(false);
            this.reset();
            setLabel("Elapsed time");
        }
        
        /**
         * Declaration that this meter does not use the boundary object of phase when making its measurements
         */
        public final boolean usesPhaseBoundary() {return false;}
        /**
         * Declaration that this meter does not use the iteratorFactory of phase when making its measurements
         */
        public final boolean usesPhaseIteratorFactory() {return false;}
        
        public Dimension getDimension() {return Dimension.TIME;}
        
        public void reset() {t0 = elapsedTime();}
        public void intervalAction(Integrator.IntervalEvent evt) {}
        public double currentValue() {return elapsedTime() - t0;}
    	public double mostRecent() {return currentValue();}
    }
    
    /** 
     * DisplayBox to present the elapsed time
     */
    public class Timer extends DisplayBox {
    
        public Timer() {this(chronoMeter());}
        public Timer(ChronoMeter meter) {
            super();
            this.setMeter(meter);
            this.setUnit(new Unit(Picosecond.UNIT));
            this.setPrecision(7);
            setSize(100,60);
        }
    }
}
    