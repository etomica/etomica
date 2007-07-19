package etomica.action.activity;

import etomica.action.Action;
import etomica.action.Activity;
import etomica.util.Arrays;

/**
 * Organizer of simulation actions to be executed in series.
 */
public class ActivityGroupSeries extends Activity implements ActivityGroup {

    public ActivityGroupSeries() {
        super();
    }
    
    /**
	 * Adds the given action to the list of actions.
	 * If action is already in list of actions to be performed, method returns
	 * without doing anything.
	 */
	public synchronized void addAction(Action newAction) {
		//FIXME this doesn't actually check that newAction isn't already in the array
		pendingActions = (Action[])Arrays.addObject(pendingActions, newAction);
		numActions++;
	}
    

    /**
     * Removes the given action from the list of actions performed by this controller.
     * If action is currently running, or is not in list of actions, method returns 
     * without doing anything.
     * @return true if the action was removed
     */
	public synchronized boolean removeAction(Action action) {
		pendingActions = (Action[]) Arrays.removeObject(pendingActions, action);
		if (pendingActions.length == numActions) {
            int oldNumCompleted = completedActions.length;
            completedActions = (Action[]) Arrays.removeObject(completedActions, action);
            return completedActions.length != oldNumCompleted;
        }
		numActions = pendingActions.length;
		return true;
	}
    
    /**
     * Marks all actions as pending.  Previously completed actions are put 
     * first in the new list of pending actions.  The Controller must be 
     * inactive or halted before calling this method.
     */
    public synchronized void reset() {
        if (currentAction != null) {
            //throw an exception?  call halt?
            return;
        }
        Action[] newPendingActions = new Action[completedActions.length+numActions];
        System.arraycopy(completedActions,0,newPendingActions,0,completedActions.length);
        System.arraycopy(pendingActions,0,newPendingActions,completedActions.length,numActions);
        numActions = newPendingActions.length;
        pendingActions = newPendingActions;
        completedActions = new Action[0];
    }
	
    /**
     * @return a list of the actions yet to be performed by this controller. 
     */
    public synchronized Action[] getPendingActions() {
        return pendingActions;
    }
    
    /**
     * @return an array containing the action currently being performed,   
     * or a zero-length array if the current action is null.
     */
    public synchronized Action[] getCurrentActions() {
        if(currentAction == null) return new Action[0];
        return new Action[] {currentAction};
     }

    /**
     * @return a list of the actions completed by this activity group.
     */
    public synchronized Action[] getCompletedActions() {
        return completedActions;
    }
    
    public synchronized Action[] getAllActions() {
    	Action[] allActions = new Action[pendingActions.length+completedActions.length
    	                                 +((currentAction==null) ? 0 : 1)];
    	System.arraycopy(completedActions,0,allActions,0,completedActions.length);
    	int i = completedActions.length;
    	if (currentAction != null) {
    		allActions[i++] = currentAction;
    	}
    	System.arraycopy(pendingActions,0,allActions,i,pendingActions.length);
    	return allActions;
    }
    
    /**
     * Causes uncompleted actions added to this group to be run in sequence.  Should not be
     * executed directly, but instead as part of the Runnable interface it is executed
     * by a thread made upon invoking the start method.
     */
    protected void run() {
        while(numActions > 0) {
            synchronized(this) {
                currentAction = pendingActions[0];
                removeAction(currentAction);
            }
            boolean exceptionThrown = false;
            currentAction.actionPerformed();

            //TODO mark this as whether completed normally
            synchronized(this) {
                completedActions = (Action[])Arrays.addObject(completedActions, currentAction);
                currentAction = null;
            }
            if(exceptionThrown || pauseAfterEachAction) {
                pauseRequested = true;
            }
             
            if(!doContinue()) {
                break;
            }
        }
    }
            
    /**
     * Requests a pause in the performance of the actions. If the current action is
     * an Activity, it is paused; if a simple Action, pause takes effect once it
     * has completed. In either case, calling thread is put in a wait state until 
     * the pause takes effect.
     */
    public synchronized void pause() {
        //        System.out.println("in AGS.pause "+isPaused()+" "+isActive());
        if (isPaused() || !isActive()) {
            return;// already paused, or not active
        }
        if (currentAction instanceof Activity) {
            //            System.out.println("pausing "+currentAction);
            ((Activity) currentAction).pause();//activity enforces pause and
                                               // has calling thread wait till
                                               // in effect
            //            System.out.println("paused "+currentAction);
        } else {//currentAction is not a pausable activity; put pause in
                // activity loop
            super.pause();
        }
    }
    
    /**
     * Removes activity group from the paused state, resuming execution where it
     * left off.
     */
    public synchronized void unPause() {
        if (currentAction instanceof Activity) {
            if (!isPaused() || !isActive()) {
                return;// not currently paused or not active
            }
            pauseRequested = false;
            ((Activity) currentAction).unPause();
        } else {
            super.unPause();
        }
    }
         
    /**
     * Request that the activity group terminate its thread as soon as possible.
     * Calling thread is caused to wait until halt is completed.
     */
    public synchronized void halt() {
        if(!isActive()) return;
        haltRequested = true;
        if(currentAction instanceof Activity) {
            ((Activity)currentAction).halt();
        }
        try {
            while(isActive()) {
                unPause();
                wait();  //make thread requesting pause wait until halt is in effect
            }
        } catch(InterruptedException e) {}
    }
    
    public synchronized boolean isPaused() {
        return super.isPaused()
                || (currentAction instanceof Activity && ((Activity) currentAction)
                        .isPaused());
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

    private static final long serialVersionUID = 1L;
    protected Action currentAction;
    protected boolean pauseAfterEachAction;
    protected int numActions;
    protected Action[] pendingActions = new Action[0];
    protected Action[] completedActions = new Action[0];
}