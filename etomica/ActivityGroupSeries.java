package etomica;

import etomica.utility.Arrays;

/**
 * Organizer of simulation actions to be executed in series.
 */
public class ActivityGroupSeries extends Activity {
    
    /**
	 * Adds the given action to the list of actions.
	 * If action is already in list of actions to be performed, method returns
	 * without doing anything.
	 */
	public synchronized void addAction(Action newAction) {
		//FIXME this doesn't actually check that newAction isn't already in the array
		actions = (Action[])Arrays.addObject(actions, newAction);
		numActions++;
	}

    /**
     * Removes the given action from the list of actions performed by this controller.
     * If action is currently running, or is not in list of actions to be performed,
     * method returns without doing anything.
     * @return true if the action was removed
     */
	public synchronized boolean removeAction(Action action) {
		actions = (Action[]) Arrays.removeObject(actions, action);
		int newNumActions = actions.length;
		if (newNumActions == numActions)
			return false;
		numActions = newNumActions;
		return true;
	}
	
    /**
     * @return a list of the actions yet to be performed by this controller. 
     */
    public synchronized Action[] pendingActions() {return actions;}
    
    /**
     * @return an array containing the action currently being performed  
     * (if there is one).
     */
    public synchronized Action[] currentActions() {return new Action[] {currentAction};}

    /**
     * @return a list of the actions completed by this controller.
     */
    public synchronized Action[] completedActions() {return completedActions;}
    
    /**
     * Causes uncompleted actions added to this controller to be run in sequence.  Should not be
     * executed directly, but instead as part of the Runnable interface it is executed
     * by a thread made upon invoking the start method.
     */
    public void run() {
    	while(numActions > 0) {
    		synchronized(this) {
    			currentAction = actions[0];
    			removeAction(currentAction);
    		}
    		boolean exceptionThrown = false;
    		try {
    			currentAction.actionPerformed();
    		}
    		catch (Exception e) {
    			//TODO write message to error stream
    			e.printStackTrace();
    			exceptionThrown = true;
    		}
    		//TODO mark this as whether completed normally
    		synchronized(this) {
    			completedActions = (Action[])Arrays.addObject(completedActions, currentAction);
    			currentAction = null;
    		}
    		if(exceptionThrown || pauseRequested || pauseAfterEachAction) doWait();
    		if(haltRequested) break;
    	}
    	synchronized(this) {
    		notifyAll();//notify any threads requesting halt and waiting for execution to complete
    	}
    }
            
    /**
     * Requests a pause in the performance of the actions. If the current action is
     * an Activity, it is paused; if a simple Action, pause takes effect once it
     * has completed. In either case, calling thread is put in a wait state until 
     * the pause takes effect.
     */
    public synchronized void pause() {
    	if(isPaused() || !isActive()) return;// already paused, or not active
        if(currentAction instanceof Activity) {
        	((Activity)currentAction).pause();//activity enforces pause and has calling thread waits till in effect
        } else {//currentAction is not a pausable activity; put pause in controller loop
	        super.pause();
        }
    }
    
    /**
     * Removes controller from the paused state, resuming execution where it left off.
     */
    public synchronized void unPause() {
    	if (!isPaused() || !isActive()) return;// not currently paused or not active
    	pauseRequested = false;
        if(currentAction instanceof Activity) {
        	((Activity)currentAction).unPause();
        } else {
    		notifyAll();
    	}
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
        unPause();
        try {
            wait();  //make thread requesting pause wait until halt is in effect
        } catch(InterruptedException e) {}
    }
    
    public boolean isPaused() {
    	return super.isPaused() || 
    			(currentAction instanceof Activity 
    					&& ((Activity)currentAction).isPaused());
    }
    
    /**
     * @return flag specifying whether activity should pause upon completing each
     * action.
     */
	public boolean isPauseAfterEachAction() {
		return pauseAfterEachAction;
	}
	/**
	 * @param pauseAfterEachAction specifies whether activity should pause upon
	 * completing each action (true), or if next action should begin immediately
	 * upon completion of current action.
	 */
	public void setPauseAfterEachAction(boolean pauseAfterEachAction) {
		this.pauseAfterEachAction = pauseAfterEachAction;
	}


    protected Action currentAction;
    protected boolean pauseAfterEachAction;
    protected boolean haltRequested;
    protected int numActions;
    protected boolean isPaused = false;
    protected Action[] actions = new Action[0];
    protected Action[] completedActions;

}//end of Controller


