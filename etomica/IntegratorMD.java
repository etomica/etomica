package etomica;

import etomica.units.*;
/**
 * Superclass of all molecular-dynamics integrators.
 * Extends the Integrator class by adding methods that 
 * set the time step.
 */

public abstract class IntegratorMD extends Integrator {
    
    public IntegratorMD(PotentialMaster potentialMaster) {
        super(potentialMaster);
        setTimeStep(Default.TIME_STEP);
        meterKE = new MeterKineticEnergy();
    }

    /**
     * Sets integration time step.
     * Updates zero-point counters used internally to manage the elapsedTime method
     */
    public void setTimeStep(double t) {
        timeStep = t;
    }
    public final double getTimeStep() {return timeStep;}
    public Dimension getTimeStepDimension() {return Dimension.TIME;}

    /**
     * reset the integrator's kinetic energy tracker
     */
    public void reset() {
        super.reset();
        meterKE.setPhase(phase);
        currentKineticEnergy = meterKE.getData();
    }

    /**
     * @return the current kinetic energy as tracked by the integrator
     */
    public double[] getKineticEnergy() {
        reutrn currentKineticEnergy;
    }
    
    /**
     * Elementary time step for the MD simulation
     */
    protected double timeStep;
    protected double[] currentKineticEnergy;
    protected MeterKineticEnergy meterKE;
    
}//end of IntegratorMD
    
