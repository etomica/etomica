package etomica.action.activity;

import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.Integrator;

public class ActivityIntegrate2 implements Activity2 {
    private final Integrator integrator;
    private final double sleepPeriod;
    private final boolean ignoreOverlap;
    public ActivityIntegrate2(Integrator integrator) {
        this(integrator, 0, false);
    }

    public ActivityIntegrate2(Integrator integrator, double sleepPeriod, boolean ignoreOverlap) {
        this.integrator = integrator;
        this.sleepPeriod = sleepPeriod;
        this.ignoreOverlap = ignoreOverlap;
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
