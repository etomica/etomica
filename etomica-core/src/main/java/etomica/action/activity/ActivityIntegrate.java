package etomica.action.activity;

import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.Integrator;

public class ActivityIntegrate implements Activity {
    private final Integrator integrator;
    private final double sleepPeriod;
    private final boolean ignoreOverlap;
    public ActivityIntegrate(Integrator integrator) {
        this(integrator, 0, false);
    }

    public ActivityIntegrate(Integrator integrator, double sleepPeriod, boolean ignoreOverlap) {
        this.integrator = integrator;
        this.sleepPeriod = sleepPeriod;
        this.ignoreOverlap = ignoreOverlap;
    }

    public Integrator getIntegrator() {
        return integrator;
    }

    @Override
    public void preAction() {
        try {
            this.integrator.reset();
        } catch (ConfigurationOverlapException e) {
            if (!ignoreOverlap) {
                throw e;
            }
        }
        integrator.resetStepCount();
    }

    @Override
    public void postAction() {

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

    @Override
    public void actionPerformed() {
        integrator.doStep();
    }

    public boolean isIgnoreOverlap() {
        return ignoreOverlap;
    }

    double getSleepPeriod() {
        return this.sleepPeriod;
    }
}
