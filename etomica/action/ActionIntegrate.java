package etomica.action;

import etomica.api.IIntegrator;
import etomica.exception.ConfigurationOverlapException;
import etomica.util.Debug;

/**
 * Action that repeatedly invokes an Integrator's doStep method.
 */
public class ActionIntegrate implements Action {

    /**
	 * Constructs activity to generate configurations with
	 * the given integrator.
	 */
	public ActionIntegrate(IIntegrator integrator) {
        this(integrator, false);
    }
    
    public ActionIntegrate(IIntegrator integrator, boolean ignoreOverlap) {
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
            integrator.initialize();
        }
        catch (ConfigurationOverlapException e) {
            if (!ignoreOverlap) {
                throw new RuntimeException(e);
            }
        }
        for (stepCount = 0; stepCount < maxSteps; stepCount++) {
            if (Debug.ON && stepCount == Debug.START) Debug.DEBUG_NOW = true;
            if (Debug.ON && stepCount == Debug.STOP) break;
            if (Debug.ON && Debug.DEBUG_NOW) System.out.println("*** integrator step "+stepCount);
            integrator.doStep();
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
	public IIntegrator getIntegrator() {
		return integrator;
	}
    
    private static final long serialVersionUID = 1L;
	private final IIntegrator integrator;
    private boolean ignoreOverlap;
	protected long maxSteps, stepCount;
}
