package etomica;

import java.util.LinkedList;

/**
 * Organizer of actions of the integrators.  Controller sets the protocol for
 * conduct of the simulation.  For example, the simulation might start and stop 
 * with the push of a button; it might run for a fixed number of cycles and exit;
 * it might run several simulations over a series of states.<p>
 * The Controller runs on its own thread, and it spawns Integrator processes
 * each on its own thread.
 */
public class Controller implements Runnable, java.io.Serializable, EtomicaElement {
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Simple controller that enables start/stop of simulation using a button");
        return info;
    }
 
   /**
    * Adds the given action to the list of actions performed by this controller.
    * If action is already in list of actions to be performed, method returns
    * without doing anything.
    */
    public synchronized void add(Action action) {
        if (!actions.contains(action)) actions.add(action);
    }

    /**
     * Removes the given action from the list of actions performed by this controller.
     * If action is currently running, or is not in list of actions to be performed,
     * method returns without doing anything.
     */
    public synchronized void remove(Action action) {
        if (action != currentAction) actions.remove(action);
    }
        
    /**
     * @return a list of the actions yet to be performed by this controller, 
     * including the action currently being performed (if there is one).
     */
    public LinkedList actions() {return actions;}
    
    /**
     * @return a list of the actions completed by this controller.
     */
    public LinkedList completedActions() {return completedActions;}
    
    /**
     * Method to start this controller's thread, causing execution of the run() method.
     * No action is performed if a thread is already running this controller (i.e., if
     * start was called before, and actions remain to be completed).
     */
    public void start() {
        if(isActive()) return;
        pauseRequested = false;
        runner = new Thread(this);
        runner.start();
    }
                    
    /**
     * Causes uncompleted actions added to this controller to be run in sequence.  Should not be
     * executed directly, but instead as part of the Runnable interface it is executed
     * by a thread made upon invoking the start method.
     */
    public void run() {
    	while(actions.size() > 0) {
    		currentAction = (Action)actions.getFirst();
    		boolean doPause = pauseAfterEachAction;
    		try {
    			currentAction.actionPerformed();
    		}
    		catch (Exception e) {
    			//TODO write message to error stream
    			e.printStackTrace();
    			doPause = true;
    		}
    		completedActions.addLast(currentAction);//TODO mark this as whether completed normally
    		actions.removeFirst();
    		currentAction = null;
    		if(doPause || pauseRequested) doWait();
    		if(haltRequested) break;
    	}
    	synchronized(this) {
    		notifyAll();
    	}
        runner = null;
    }
            
    /**
     * Returns true if a thread exists performing this controller's actions.
     */
    public boolean isActive() {
    	return runner != null;
    }
    
    /**
     * Method to put controller in a condition of being paused.
     */
    private synchronized void doWait() {
        notifyAll(); //release any threads waiting for pause to take effect
        try {
            wait(); //put in paused state
        } catch(InterruptedException e) {}
        pauseRequested = false;
    }
    
    /**
     * Requests a pause in the performance of the actions. If the current action is
     * an Activity, it is paused; if a simple Action, pause takes effect once it
     * has completed. In either case, calling thread is put in a wait state until 
     * the pause takes effect.
     */
    public synchronized void pause() {
        if(!isActive() || currentAction == null) return;//nothing going on, no need to pause
        
        if(currentAction instanceof Activity) {
        	((Activity)currentAction).pause();//activity enforces pause and has calling thread waits till in effect
        } else {//currentAction is not a pausable activity; put pause in controller loop
        	pauseRequested = true;
	        try {
	            wait();  //make thread requesting pause wait until pause is in effect
	        } catch(InterruptedException e) {}
        }
    }
    
    /**
     * Removes controller from the paused state, resuming execution where it left off.
     */
    public synchronized void unPause() {
        if(currentAction != null && currentAction instanceof Activity) {
        	((Activity)currentAction).unPause();
        } else {
    		notifyAll();
    	}
    	pauseRequested = false;

    }
    
    /**
     * Queries whether the integrator is in a state of being paused.  This may
     * occur independent of whether the integrator is running or not.  If paused
     * but not running, then pause will take effect upon start.
     */
    public boolean isPaused() {
    	if(currentAction instanceof Activity) {
    		return ((Activity)currentAction).isPaused();
    	} else {
    		return pauseRequested;
    	}
    }
     
    /**
     * @return flag specifying whether controller should pause upon completing each
     * action.
     */
	public boolean isPauseAfterEachActivity() {
		return pauseAfterEachAction;
	}
	/**
	 * @param pauseAfterEachAction specifies whether controller should pause upon
	 * completing each action (true), or if next action should begin immediately
	 * upon completion of current action.
	 */
	public void setPauseAfterEachActivity(boolean pauseAfterEachAction) {
		this.pauseAfterEachAction = pauseAfterEachAction;
	}

    /**
     * Request that the integrator terminate its thread on the next integration step.
     * Does not cause calling thread to wait until this is completed, so it would
     * be prudent to have the calling thread join() to suspend it until the halt
     * is in effect.
     */
    public synchronized void halt() {
        if(!isActive()) return;
        haltRequested = true;
        if(currentAction instanceof Activity) ((Activity)currentAction).halt();
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
    
    private Action currentAction;
    private boolean pauseAfterEachAction;
    private boolean pauseRequested;
    private boolean haltRequested;
    private final LinkedList actions = new LinkedList();
    private final LinkedList completedActions = new LinkedList();
    private transient Thread runner;      
    private SimulationEventManager eventManager = new SimulationEventManager();

}//end of Controller


