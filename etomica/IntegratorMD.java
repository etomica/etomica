package etomica;

import etomica.units.*;
/**
 * Superclass of all molecular-dynamics integrators.
 * Extends the Integrator class by adding methods that 
 * set the time step.
 */

//TODO consider making an interface

public abstract class IntegratorMD extends Integrator {
    
    public IntegratorMD(PotentialMaster potentialMaster) {
        super(potentialMaster);
        setTimeStep(Default.TIME_STEP);
    }
    
    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return false;
        return true;
    }

    /**
     * Sets integration time step.
     * Updates zero-point counters used internally to manage the elapsedTime method
     */
    public void setTimeStep(double t) {
        timeStep = t;
    }
    public final double getTimeStep() {return timeStep;}
    public final double timeStep() {return timeStep;}
    public Dimension getTimeStepDimension() {return Dimension.TIME;}
    
    /**
     * Elementary time step for the MD simulation
     */
    protected double timeStep;
    
}//end of IntegratorMD
    