package etomica.action.activity;

import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.Integrator;
import etomica.util.Debug;

public class ActivityIntegrate extends Activity {

    private final Integrator integrator;

    private long maxSteps;
    private boolean ignoreOverlap;
    private long currentStep = 0;

    public ActivityIntegrate(Integrator integrator, long maxSteps, boolean ignoreOverlap) {
        this.integrator = integrator;
        this.maxSteps = maxSteps;
        this.ignoreOverlap = ignoreOverlap;
    }

    public ActivityIntegrate(Integrator integrator, long maxSteps) {
        this(integrator, maxSteps, false);
    }

    public ActivityIntegrate(Integrator integrator, boolean ignoreOverlap) {
        this(integrator, Long.MAX_VALUE, ignoreOverlap);
    }

    public ActivityIntegrate(Integrator integrator) {
        this(integrator, Long.MAX_VALUE, false);
    }

    public void doAction() {
        if (Debug.ON) {
            if (currentStep == Debug.START) { Debug.DEBUG_NOW = true; }
            if (Debug.DEBUG_NOW && Debug.LEVEL > 0) {
                System.out.println("*** integrator step " + currentStep);
            }
            Debug.stepCount = currentStep;
        }
        this.integrator.doStep();
    }


    @Override
    public void runActivity(Controller.ControllerHandle handle) {
        try {
            this.integrator.reset();
        } catch (ConfigurationOverlapException e) {
            if (!ignoreOverlap) {
                throw e;
            }
        }
        integrator.resetStepCount();

        for (currentStep = 0; currentStep < this.maxSteps; currentStep++) {
            handle.yield(this::doAction);
            if (Debug.ON && currentStep == Debug.STOP) { break; }
        }

    }

    @Override
    public void restart() {
        if (integrator.getStepCount() > 0) {
            integrator.resetStepCount();
        }
        if (integrator.isInitialized()) {
            try {
                integrator.reset();
            } catch (ConfigurationOverlapException e) {
                if (!ignoreOverlap) {
                    throw e;
                }
            }
        }
    }

    public Integrator getIntegrator() {
        return integrator;
    }

    public long getMaxSteps() {
        return maxSteps;
    }

    public void setMaxSteps(long maxSteps) {
        this.maxSteps = maxSteps;
    }

    public boolean isIgnoreOverlap() {
        return ignoreOverlap;
    }

    public void setIgnoreOverlap(boolean ignoreOverlap) {
        this.ignoreOverlap = ignoreOverlap;
    }

}
