/*
 * History
 * Created on Oct 26, 2004 by kofke
 */
package etomica;

import etomica.Integrator.IntervalEvent;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ActivityIntegrate extends Activity {

	/**
	 * 
	 */
	public ActivityIntegrate(Integrator integrator) {
		super();
		this.integrator = integrator;
		interval = 1;
		doSleep = Default.DO_SLEEP;
		sleepPeriod = 10;
	}

    /**
     * Main loop for conduct of integration.  Repeatedly calls doStep() method,
     * while checking for halt/pause/reset requests, firing regular interval events,
     * and entering a brief sleep state if so indicated by doSleep flag.  Integration
     * loop continues until number of steps equals maxSteps field.  This method
     * is not the normal way to run the integrator, as it runs on the calling thread.
     * Instead, start should be used to set the integrator to run on its own new thread.
     */
	public void run() {
        integrator.fireIntervalEvent(new IntervalEvent(integrator, IntervalEvent.START));
	    integrator.initialize();
		IntervalEvent intervalEvent = new IntervalEvent(integrator, IntervalEvent.INTERVAL);
        int stepCount = 0;
        int iieCount = interval;//changed from "interval + 1"
        while(stepCount < maxSteps) {
        	if (Debug.ON && stepCount == Debug.START) Debug.DEBUG_NOW = true;
        	if (Debug.ON && stepCount == Debug.STOP) halt();
        	if (Debug.ON && Debug.DEBUG_NOW) System.out.println("*** integrator step "+stepCount);
        	if (!doContinue()) break;
        	if (integrator.resetRequested()) {integrator.doReset();}
            integrator.doStep();
            if(--iieCount == 0) {
                integrator.fireIntervalEvent(intervalEvent);
                iieCount = interval;
            }
            if(doSleep) {
                try { Thread.sleep(sleepPeriod); }
                catch (InterruptedException e) { }
            }
            stepCount++;
        }//end of while loop
        integrator.fireIntervalEvent(new IntervalEvent(integrator, IntervalEvent.DONE));
	}

	public boolean isDoSleep() {
		return doSleep;
	}
	public void setDoSleep(boolean doSleep) {
		this.doSleep = doSleep;
	}
	public int getInterval() {
		return interval;
	}
	public void setInterval(int interval) {
		this.interval = interval;
	}
	public int getSleepPeriod() {
		return sleepPeriod;
	}
	public void setSleepPeriod(int sleepPeriod) {
		this.sleepPeriod = sleepPeriod;
	}
    /**
     * Accessor method for the number of integration steps to be
     * performed by this integrator after it is started.
     */
	public int getMaxSteps() {
		return maxSteps;
	}

	/**
     * Mutator method for the number of integration steps to be
     * performed by this integrator after it is started.  Can
     * be changed while integrator is running; if set to a value
     * less than number of steps already executed, integration will end.
     */
	public void setMaxSteps(int maxSteps) {
		this.maxSteps = maxSteps;
	}

	private final Integrator integrator;
	protected int interval;
	private boolean doSleep;
	private int sleepPeriod;
	protected int maxSteps = Integer.MAX_VALUE;
}
