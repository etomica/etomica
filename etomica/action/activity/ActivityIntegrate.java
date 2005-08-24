/*
 * History
 * Created on Oct 26, 2004 by kofke
 */
package etomica.action.activity;

import etomica.action.Activity;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.util.Debug;
import etomica.util.Default;

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
        setInterval(1);
        doSleep = Default.DO_SLEEP;
        sleepPeriod = 10;
        setLabel("ActivityIntegrate");
	}
    
    protected ActivityIntegrate(ActivityIntegrate activity) {
        super(activity);
        integrator = activity.integrator;
        setDoSleep(activity.doSleep);
        setInterval(activity.interval);
        setMaxSteps(activity.maxSteps);
        setSleepPeriod(activity.sleepPeriod);
    }
    
    public Activity makeCopy() {
        return new ActivityIntegrate(this);
    }

    /**
     * Main loop for conduct of integration.  Repeatedly calls doStep() method,
     * while checking for halt/pause/reset requests, firing regular interval events,
     * and entering a brief sleep state if so indicated by doSleep flag.  Integration
     * loop continues until number of steps equals maxSteps field.  This method should
     * not be called directly, but instead is called by the instance's actionPerformed method.
     */
	public void run() {
        integrator.fireNonintervalEvent(new IntegratorNonintervalEvent(integrator, IntegratorNonintervalEvent.START));
	    integrator.initialize();
        long stepCount = 0;
        int iieCount = interval;//changed from "interval + 1"
        while(stepCount < maxSteps) {
        	if (Debug.ON && stepCount == Debug.START) Debug.DEBUG_NOW = true;
        	if (Debug.ON && stepCount == Debug.STOP) halt();
        	if (Debug.ON && Debug.DEBUG_NOW) System.out.println("*** integrator step "+stepCount);
        	if (!doContinue()) break;
        	if (resetRequested) {integrator.reset();}
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
        integrator.fireNonintervalEvent(new IntegratorNonintervalEvent(integrator, IntegratorIntervalEvent.DONE));
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
		intervalEvent = new IntegratorIntervalEvent(integrator, interval);
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
		if(maxSteps < 0) maxSteps = 0;
		this.maxSteps = maxSteps;
	}
	
	/**
	 * Requests that the integrator reset itself; if integrator is not running,
	 * then reset will be performed immediately. The actual action an integrator
	 * takes to do this differs with the type of integrator. The reset is not
	 * performed until the completion of the current integration step, or until
	 * the integrator is unpaused if currently in a paused state.
	 */
	public void reset() {
		resetRequested = true;
	}

	public boolean resetRequested() {
		return resetRequested;
	}


	/**
	 * @return Returns the integrator.
	 */
	public Integrator getIntegrator() {
		return integrator;
	}
	
	private final Integrator integrator;
	protected int interval;
	private boolean resetRequested;
	private boolean doSleep;
	private int sleepPeriod;
	protected long maxSteps = Long.MAX_VALUE;
	IntegratorIntervalEvent intervalEvent;

}
