/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.integrator.Integrator;
import etomica.exception.ConfigurationOverlapException;
import etomica.util.Debug;

/**
 * Action that repeatedly invokes an Integrator's doStep method.
 */
public class ActionIntegrate implements IAction {

    /**
	 * Constructs activity to generate configurations with
	 * the given integrator.
	 */
	public ActionIntegrate(Integrator integrator) {
        this(integrator, false);
    }
    
    public ActionIntegrate(Integrator integrator, boolean ignoreOverlap) {
        super();
        this.integrator = integrator;
        this.ignoreOverlap = ignoreOverlap;
        setMaxSteps(Integer.MAX_VALUE);
	}
    
    /**
     * Main loop for conduct of integration.  Repeatedly calls doStep() method,
     * while checking for halt/pause/reset requests, firing regular interval events,
     * and entering a brief sleep state if so indicated by doSleep flag.  Integration
     * loop continues until number of steps equals maxSteps field.  This method should
     * not be called directly, but instead is called by the instance's actionPerformed method.
     */
    public void actionPerformed() {
        try {
            integrator.reset();
        }
        catch (ConfigurationOverlapException e) {
            if (!ignoreOverlap) {
                throw e;
            }
        }
        integrator.resetStepCount();
        for (stepCount = 0; stepCount < maxSteps; stepCount++) {
            long t0 = System.nanoTime();
            if (Debug.ON) {
                if (stepCount == Debug.START) Debug.DEBUG_NOW = true;
                if (stepCount == Debug.STOP) break;
                if (Debug.DEBUG_NOW && Debug.LEVEL > 0) System.out.println("*** integrator step "+stepCount);
                Debug.stepCount = stepCount;
            }
            integrator.doStep();
            long t1 = System.nanoTime();
            System.out.println("Step " + stepCount + " in " + (t1 - t0) / 1_000_000 + " ms");
        }
	}

    /**
     * Accessor method for the number of doStep calls to be
     * performed by this integrator after it is started.
     */
	public long getMaxSteps() {
		return maxSteps;
	}

	/**
     * Mutator method for the number of doStep steps to be
     * performed by this integrator after it is started.  Can
     * be changed while activity is running; if set to a value
     * less than number of steps already executed, integration will end.
     */
	public void setMaxSteps(long maxSteps) {
		if(maxSteps < 0) throw new IllegalArgumentException("steps must not be negative");
		this.maxSteps = maxSteps;
	}
    
    public long getCurrentStep() {
        return stepCount;
    }
	
	/**
	 * @return Returns the integrator.
	 */
	public Integrator getIntegrator() {
		return integrator;
	}
    
    private static final long serialVersionUID = 1L;
	private final Integrator integrator;
    private boolean ignoreOverlap;
	protected long maxSteps, stepCount;
}
