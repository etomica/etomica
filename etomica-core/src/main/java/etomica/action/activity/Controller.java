/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action.activity;

import java.util.HashMap;

import etomica.action.Activity;
import etomica.action.IAction;
import etomica.util.Arrays;
import etomica.util.EnumeratedType;
import etomica.util.IEventManager;

/**
 * Organizer and executor of actions performed by the simulation. The Controller
 * holds a series of Actions which it performs in series, using a thread that it
 * launches when start() is invoked. A simple Action will run on the
 * Controller's thread. An Activity is a subclass of Action that runs on its own
 * thread, and can be paused or terminated (halted) before it completes
 * naturally.
 * <p>
 * The Controller thread will idle while an Activity is running. While idling,
 * the Controller will respond to "urgent actions" that another thread
 * (typically a gui) will ask it to perform. In response, the Controller will
 * pause the current Activity (or complete the Action if not an Activity),
 * perform the urgent action, and continue where it left off. This capability is
 * needed for a user to interact with a simulation that is in progress, without
 * colliding with an integrator running on this Controller.
 * <p>
 * The queue for pending Actions may be altered (given additions or deletions)
 * while Actions are being executed. As Actions are completed they are moved
 * into a list of completed Actions. When all Actions are done, the Controller
 * exits, but it may be restarted if Actions are added and start() is invoked
 * again.
 * <p>
 * Each Simulation holds a single, final instance of a Controller.
 * 
 * @author Andrew Schultz and David Kofke
 * 
 * @see IAction
 * @see Activity
 */
public class Controller extends ActivityGroupSeries {
    
    public Controller() {
        actionStatusMap = new HashMap<IAction, ActionStatus>();
        actionExceptionMap = new HashMap<IAction, Throwable>();
        waitObject = new WaitObject();
        urgentWaitObject = new UrgentWaitObject();
        eventManager = new ControllerEventManager();
    }
    
    public synchronized void addAction(IAction newAction) {
        super.addAction(newAction);
        actionStatusMap.put(newAction,ActionStatus.PENDING);
        // don't need to add to the exception map because HashMap will return
        // null if we don't add it.
    }
    
    public synchronized boolean removeAction(IAction oldAction) {
        if (super.removeAction(oldAction)) {
            actionStatusMap.remove(oldAction);
            actionExceptionMap.remove(oldAction);
            return true;
        }
        return false;
    }

    /**
     * Marks all actions as pending.  Previously completed actions are put
     * first in the new list of pending actions.  The Controller must be
     * inactive or halted before calling this method.
     */
    public synchronized void reset() {
        for (int i=0; i<completedActions.length; i++) {
            actionStatusMap.remove(completedActions[i]);
            actionStatusMap.put(completedActions[i],ActionStatus.PENDING);
            actionExceptionMap.remove(completedActions[i]);
        }
        super.reset();
        eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.RESET, currentAction));
    }
    
    /**
     * Returns the status of an action held by the controller.  Returns
     * null for actions not held by the controller.
     */
    public synchronized ActionStatus getActionStatus(IAction action) {
        return actionStatusMap.get(action);
    }
    
    /**
     * Returns the exception thrown by an action held by the controller.  
     * Returns null if the given action did not throw an exception or is not 
     * held by the controller.
     */
    public synchronized Throwable getException(IAction action) {
        return actionExceptionMap.get(action);
    }
    
    /**
     * Causes uncompleted actions added to this group to be run in sequence.
     * Should not be executed directly, but instead it is executed by a thread
     * made upon invoking the start method.
     */
    protected void run() {
        Thread.currentThread().setName("Controller thread");
        eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.START));
        while(numActions > 0) {
            synchronized(this) {
                currentAction = pendingActions[0];
                //removeAction will remove the action from the status map
                removeAction(currentAction);
                //put it back
                actionStatusMap.put(currentAction,ActionStatus.CURRENT);
            }
            eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.START_ACTION, currentAction));
            if(currentAction instanceof Activity) {
                waitObject.currentActionDone = false;
                waitObject.actionException = null;
                //define a thread to run the activity
                Thread localRunner = new Thread(new Runnable() {
                    public void run() {
                        Thread.currentThread().setName(currentAction+" thread");
                        try {
                            currentAction.actionPerformed();
                        }
                        catch (Throwable e) {
                            synchronized(waitObject) {
                                waitObject.actionException = e;
                            }
                        }
                        synchronized(waitObject) {
                            waitObject.currentActionDone = true;
                            waitObject.notifyAll();
                        }
                    }
                });
                //start the activity's thread, and wait for it to exit
                localRunner.start();
                while(!waitObject.currentActionDone) {
                    try {
                        synchronized(waitObject) {
                            //put group thread in wait state while current action runs
                            //need to check again that current action isn't done, to ensure that 
                            //current-action thread didn't finish between while and here
                            while(!waitObject.currentActionDone) {
                                waitObject.wait();
                            }
                        }
                    } catch(InterruptedException e) {}

                }
            } else {//currentAction is not an Activity; run on group's thread
                currentAction.actionPerformed();
            }
            
            synchronized(this) {
                //update the action's status
                actionStatusMap.remove(currentAction);
                if (waitObject.actionException != null) {
                    actionStatusMap.put(currentAction,ActionStatus.FAILED);
                    actionExceptionMap.put(currentAction,waitObject.actionException);
                }
                else if (haltRequested) {
                    actionStatusMap.put(currentAction,ActionStatus.STOPPED);
                }
                else {
                    actionStatusMap.put(currentAction,ActionStatus.COMPLETED);
                }
                completedActions = (IAction[])Arrays.addObject(completedActions, currentAction);
                eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.END_ACTION, currentAction));

                currentAction = null;
            }

            if (waitObject.actionException != null) {
                eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.NO_MORE_ACTIONS));
                isActive = false;
                throw new RuntimeException("action failed", waitObject.actionException);
            }
            
            if(pauseAfterEachAction) {
                pauseRequested = true;
            }
            
            if(haltRequested) {
                eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.HALTED));
                haltRequested = false;
            }
            
            if(!doContinue()) {
                break;
            }

        }//end while(numActions > 0)
        isActive = false;
        synchronized(this) {
            notifyAll();//notify any threads requesting halt and waiting for execution to complete
        }

        eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.NO_MORE_ACTIONS));
    }

    /**
     * Pauses current activity, executes given action, then resumes current
     * activity. If current action is already paused, it is not resumed by this
     * method.
     * 
     * @param action
     *            Action to be performed right away; cannot be an Activity.
     */
    public synchronized void doActionNow(final IAction action) {
        if(action instanceof Activity) throw new IllegalArgumentException("can't send Activity here; sorry");

        // use extra synchronization to block calls to doActionNow while we're 
        // in pause and unPause
        synchronized (urgentWaitObject) {
            if(isActive()) {
                //If the controller is active, pause it
                boolean wasPaused = isPaused();
                if (!wasPaused) {
                    pause();
                }
                
                // catch any exception this throws and then rethrow it after
                // we've unpaused ourselves
                RuntimeException thrownException = null;
                try {
                    doUrgentAction(action);
                }
                catch (RuntimeException e) {
                    thrownException = e;
                }

                //If we had to pause the controller, unpause it now
                if (!wasPaused) {
                    unPause();
                }
                
                if (thrownException != null) {
                    throw thrownException;
                }
            }
            else {
                doUrgentAction(action);
            }
        }
    }
    
    private void doUrgentAction(IAction urgentAction) {
        eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.START_URGENT_ACTION, urgentAction));
        
        urgentAction.actionPerformed();

        eventManager.fireEvent(new ControllerEvent(this, ControllerEvent.END_URGENT_ACTION, urgentAction));
    }
    
    /**
     * Flag indicating whether the current action should be repeated indefinitely before
     * moving on to the next action.  Setting this to true also sets pauseAfterEachAction 
     * to true (setting false then has no effect on pauseAfterEachAction).  In repeating
     * an action, when an action is completed it is placed again at the head of the pending 
     * action list (for an Action the same instance is used; for an Activity a copy is made).
     * <p>
     * Default is false.
     */
    public synchronized void setRepeatCurrentAction(boolean repeatCurrentAction) {
        this.repeatCurrentAction = repeatCurrentAction;
        if(repeatCurrentAction) {
            setPauseAfterEachAction(true);
        }
    }
    
    /**
     * Returns value of repeatCurrentAction flag.
     */
    public boolean isRepeatCurrentAction() {
        return repeatCurrentAction;
    }

    public String toString() {
        return "Controller";
    }

    /**
     * Returns the event manager used by the controller to notify listeners of
     * individual actions events as well as events related to the controller
     * itself.
     * @see ControllerEvent
     */
    public IEventManager getEventManager() {
        return eventManager;
    }
    
    private final IEventManager eventManager;

    private boolean repeatCurrentAction = false;

    protected final UrgentWaitObject urgentWaitObject;

    protected static class UrgentWaitObject implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
    }
    
    protected final WaitObject waitObject;
    protected static class WaitObject implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        boolean currentActionDone;
        Throwable actionException;
    };

    private static final long serialVersionUID = 1L;
    protected final HashMap<IAction,ActionStatus> actionStatusMap;
    protected final HashMap<IAction,Throwable> actionExceptionMap;
    
    /**
     * Enumerated type describing the status of an action.
     */
    public static class ActionStatus extends EnumeratedType {
        protected ActionStatus(String label) {
            super(label);
        }       

        public static final ActionStatus PENDING = new ActionStatus("Pending");
        public static final ActionStatus CURRENT = new ActionStatus("Current");
        public static final ActionStatus COMPLETED = new ActionStatus("Completed");
        public static final ActionStatus STOPPED = new ActionStatus("Stopped");
        public static final ActionStatus FAILED = new ActionStatus("Failed");

        public ActionStatus[] choices() { 
            return new ActionStatus[] {PENDING,CURRENT,COMPLETED,STOPPED,FAILED};
        }
        
        /**
         * Required to guarantee singleton when deserializing.
         * @return the singleton INSTANCE
         */
        private Object readResolve() {
            ActionStatus[] choices = choices();
            for (int i=0; i<choices.length; i++) {
                if (this.toString().equals(choices[i].toString())) {
                    return choices[i];
                }
            }
            throw new RuntimeException("unknown action status: "+this);
        }

        private static final long serialVersionUID = 1L;
    }

}
