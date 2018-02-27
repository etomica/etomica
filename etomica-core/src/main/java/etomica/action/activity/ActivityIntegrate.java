/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action.activity;

import etomica.action.Activity;
import etomica.integrator.Integrator;
import etomica.exception.ConfigurationOverlapException;
import etomica.util.Debug;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Activity that repeatedly invokes an Integrator's doStep method.
 */
public class ActivityIntegrate extends Activity {
    private final Logger LOG = LogManager.getLogger(this.getClass());

    /**
	 * Constructs activity to generate configurations with
	 * the given integrator (which is final).  Defaults include
	 * interval = 1, doSleep given by Default class, and sleepPeriod = 10.
	 */
	public ActivityIntegrate(Integrator integrator) {
        this(integrator,0,false);
    }
    
    public ActivityIntegrate(Integrator integrator, int sleepPeriod, boolean ignoreOverlap) {
        super();
        this.integrator = integrator;
        this.ignoreOverlap = ignoreOverlap;
        this.sleepPeriod = sleepPeriod;
        setMaxSteps(Long.MAX_VALUE);
	}
    
    /**
     * Main loop for conduct of integration.  Repeatedly calls doStep() method,
     * while checking for halt/pause/reset requests, firing regular interval events,
     * and entering a brief sleep state if so indicated by doSleep flag.  Integration
     * loop continues until number of steps equals maxSteps field.  This method should
     * not be called directly, but instead is called by the instance's actionPerformed method.
     */
    protected void run() {
        try {
            integrator.reset();
        }
        catch (ConfigurationOverlapException e) {
            if (!ignoreOverlap) {
                throw e;
            }
        }
        integrator.resetStepCount();
        long startTime = 0, endTime;
        for (stepCount = 0; stepCount < maxSteps; stepCount++) {
            if (LOG.isDebugEnabled()) {
                startTime = System.nanoTime();
                LOG.debug("Beginning step {}", stepCount);
            }
            if (Debug.ON) {
                if (stepCount == Debug.START) Debug.DEBUG_NOW = true;
                if (stepCount == Debug.STOP) break;
                if (Debug.DEBUG_NOW && Debug.LEVEL > 0) System.out.println("*** integrator step "+stepCount);
                Debug.stepCount = stepCount;
            }
            if (!doContinue()) break;
            integrator.doStep();
            if(sleepPeriod > 0) {
                try { Thread.sleep(sleepPeriod); }
                catch (InterruptedException e) { }
            }

            if (LOG.isDebugEnabled()) {
                endTime = System.nanoTime();
                LOG.debug("Step {} done in {}ms", stepCount, (endTime - startTime) / 1_000_000);
            }
        }
	}

	/**
	 * Amount of time that thread is kept in sleep state after
	 * each doStep done on integrator.  If doSleep is false, this no sleep
	 * is performed and this parameter has no effect.
	 * 
	 * @return sleep period, in milliseconds.
	 */
	public int getSleepPeriod() {
		return sleepPeriod;
	}
	
	/**
	 * Sets amount of time that thread is kept in sleep state after
	 * each doStep done on integrator.  Default value is 0.
	 */

	public void setSleepPeriod(int sleepPeriod) {
		this.sleepPeriod = sleepPeriod;
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

    @Override
    public String toString() {
        return "integration";
    }

    private static final long serialVersionUID = 1L;
	private final Integrator integrator;
    private boolean ignoreOverlap;
	private int sleepPeriod;
	protected long maxSteps, stepCount;
}
