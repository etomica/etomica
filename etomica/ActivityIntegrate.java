/*
 * History
 * Created on Oct 26, 2004 by kofke
 */
package etomica;

import etomica.Integrator.IntervalEvent;

/**
 * Activity that repeatedly invokes an Integrator's doStep method.
 */
public class ActivityIntegrate extends Activity {

	/**
	 * Constructs activity to generate configurations with
	 * the given integrator (which is final).  Defaults include
	 * interval = 1, doSleep given by Default class, and sleepPeriod = 10.
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
     * loop continues until number of steps equals maxSteps field.  This method should
     * not be called directly, but instead is called by the instance's actionPerformed method.
     */
	public void run() {
        integrator.fireIntervalEvent(new IntervalEvent(integrator, IntervalEvent.START));
	    integrator.initialize();
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

	/**
	 * @return flag specifying whether Activity puts thread to sleep briefly after each
	 * call to integrator's doStep.  This might be desired when performing interactive
	 * simulations, so the graphical interface is given time to update and respond to
	 * user input.
	 */
	public boolean isDoSleep() {
		return doSleep;
	}
	/**
	 * Sets flag specifying whether Activity puts thread to sleep briefly after each
	 * call to integrator's doStep. This might be desired when performing interactive
	 * simulations, so the graphical interface is given time to update and respond to
	 * user input.  For any type of production calculation, this should be set to false,
	 * as it seriously hampers performance.
	 */
	public void setDoSleep(boolean doSleep) {
		this.doSleep = doSleep;
	}
	
	/**
	 * @return value of interval (number of doStep calls) between
	 * firing of interval events by integrator.
	 */
	public int getInterval() {
		return interval;
	}
	
	/**
	 * Sets value of interval between successive firing of integrator interval events.
	 * @param interval
	 */
	public void setInterval(int interval) {
		this.interval = interval;
		intervalEvent = new IntervalEvent(integrator, interval);
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
	 * each doStep done on integrator.  If doSleep is false, this no sleep
	 * is performed and this parameter has no effect.  Default value is 10.
	 */

	public void setSleepPeriod(int sleepPeriod) {
		this.sleepPeriod = sleepPeriod;
	}
    /**
     * Accessor method for the number of doStep calls to be
     * performed by this integrator after it is started.
     */
	public int getMaxSteps() {
		return maxSteps;
	}

	/**
     * Mutator method for the number of doStep steps to be
     * performed by this integrator after it is started.  Can
     * be changed while activity is running; if set to a value
     * less than number of steps already executed, integration will end.
     */
	public void setMaxSteps(int maxSteps) {
		if(maxSteps < 0) maxSteps = 0;
		this.maxSteps = maxSteps;
	}

	private final Integrator integrator;
	protected int interval;
	private boolean doSleep;
	private int sleepPeriod;
	protected int maxSteps = Integer.MAX_VALUE;
	IntervalEvent intervalEvent;

}
