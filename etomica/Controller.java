package etomica;

//Java2 imports
import java.util.LinkedList;

//import etomica.utility.java2.LinkedList;
//import etomica.utility.java2.Iterator;

/**
 * Organizer of actions of the integrators.  Controller sets the protocol for
 * conduct of the simulation.  For example, the simulation might start and stop 
 * with the push of a button; it might run for a fixed number of cycles and exit;
 * it might run several simulations over a series of states.<p>
 * The Controller runs on its own thread, and it spawns Integrator processes
 * each on its own thread.
 */
public class Controller implements Runnable, java.io.Serializable, EtomicaElement {

  /**
   * List of integrators managed by the controller
   */
    private final LinkedList actions = new LinkedList();
    private final LinkedList completedActions = new LinkedList();
  /**
   * Thread used to run the controller
   */
    protected transient Thread runner;
    
    private SimulationEventManager eventManager = new SimulationEventManager();

    /**
     * Creates a new controller.
     *
     */
    public Controller() {
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Simple controller that enables start/stop of simulation using a button");
        return info;
    }
 
   /**
    * Adds an integrator to the list of integrators managed by this controller.
    * Starts the integrator immediately if autoStart is true
    */
    public synchronized void add(Action i) {
        if (!actions.contains(i)) actions.add(i);
    }

    public synchronized void remove(Action i) {
        if (i != currentAction) actions.remove(i);
    }
        
    /**
     * @return a list of the integrators managed by this controller
     */
    public LinkedList actions() {return actions;}
    
    /**
     * Accessor method for the autoStart flag.
     * autoStart flags for whether the integrators start as soon as added to the controller.
     * Default is false
     */
    public void setAutoStart(boolean b) {autoStart = b;}
    /**
     * Accessor method for the autoStart flag.
     * autoStart flags for whether the integrators start as soon as added to the controller.
     * Default is false
     */
    public boolean getAutoStart() {return autoStart;}
    

    /**
     * Method to start this controller's thread, causing execution of the run() method
     */
    public void start() {
        if(runner != null) return;
        paused = false;
        runner = new Thread(this);
        runner.start();
    }
                    
    /**
     * Activity performed by this controller once its start method is called.
     * This controller simply initiates each integrator on its own thread, which
     * is initiated by calling its start method.  Controller subclasses would override
     * this method to perform some other function.<p>
     * An integrator can be run on the controller's thread by instead calling
     * the run() method of the integrator.
     */
    public void run() {
    	while(actions.size() > 0) {
    		currentAction = (Action)actions.getFirst();
    		boolean doPause = pauseAfterEachAction;
    		try {
    			currentAction.actionPerformed();
    		}
    		catch (Exception e) {
    			doPause = true;
    		}
    		completedActions.addLast(currentAction);
    		actions.removeFirst();
    		currentAction = null;
    		if(haltRequested) break;
    		if(doPause) pause();
    	}
    	synchronized(this) {
    		notifyAll();
    	}
        runner = null;
    }
            
    /**
     * Returns true if any integrator governed by this controller is active.
     * Returns false if none are active.
     */
    public boolean isActive() {
    	return runner != null;
    }
    
    /**
     * Method to put controller in a condition of being paused.
     */
    private synchronized void doWait() {
//        isPaused = true;
		//System.out.println("pausing");
        notifyAll(); //release any threads waiting for pause to take effect
        try {
            wait(); //put in paused state
        } catch(InterruptedException e) {}
//        isPaused = false;
		//System.out.println("done pausing");
    }
    
    //suspend and resume functions
    /**
     * Requests that the integrator pause its execution.  The actual suspension
     * of execution occurs only after completion of the current integration step.
     * The calling thread is put in a wait state until the pause takes effect.
     */
    public synchronized void pause() {
        if(isActive() && !pauseRequested/*!isPaused*/) {
            pauseRequested = true;
            try {
                wait();  //make thread requesting pause wait until pause is in effect
            } catch(InterruptedException e) {}
        }
    }
    /**
     * Removes the integrator from the paused state, resuming execution where it left off.
     */
    public synchronized void unPause() {pauseRequested = false; notifyAll();}
    /**
     * Queries whether the integrator is in a state of being paused.  This may
     * occur independent of whether the integrator is running or not.  If paused
     * but not running, then pause will take effect upon start.
     */
    public boolean isPaused() {return pauseRequested;}//isPaused;}
     
	public boolean isPauseAfterEachActivity() {
		return pauseAfterEachAction;
	}
	public void setPauseAfterEachActivity(boolean pauseAfterEachAction) {
		this.pauseAfterEachAction = pauseAfterEachAction;
	}

    //stop function
    //consider having calling thread here join() or wait() for halt to take effect
    /**
     * Request that the integrator terminate its thread on the next integration step.
     * Does not cause calling thread to wait until this is completed, so it would
     * be prudent to have the calling thread join() to suspend it until the halt
     * is in effect.
     */
    public synchronized void halt() {
        if(isActive()) haltRequested = true;
        if(currentAction != null && currentAction instanceof Activity) ((Activity)currentAction).halt();
        if(pauseRequested) unPause();
        try {
            wait();  //make thread requesting pause wait until halt is in effect
        } catch(InterruptedException e) {}
    }
        
    
    public synchronized void addListener(ControllerListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(ControllerListener listener) {
        eventManager.removeListener(listener);
    }

    protected void fireEvent(ControllerEvent event) {
        eventManager.fireEvent(event);
    }    
    
    private boolean initialized = false;
    private boolean autoStart = false;
    private int maxSteps;
    private boolean paused = true;
    private Action currentAction;
    private boolean pauseAfterEachAction;
    private boolean pauseRequested;
    private boolean haltRequested;

}//end of Controller


