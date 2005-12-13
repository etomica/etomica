package etomica.action.activity;

import java.util.HashMap;

import etomica.action.Action;
import etomica.action.Activity;
import etomica.util.Arrays;
import etomica.util.EnumeratedType;

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
 * @see Action
 * @see Activity
 */
public class Controller extends ActivityGroupSeries implements java.io.Serializable {
    
    public Controller() {
        actionStatusMap = new HashMap();
        actionExceptionMap = new HashMap();
        waitObject = new WaitObject();
        eventManager = new ControllerEventManager();
    }
    
    public synchronized void addAction(Action newAction) {
        super.addAction(newAction);
        actionStatusMap.put(newAction,PENDING);
        // don't need to add to the exception map because HashMap will return
        // null if we don't add it.
    }
    
    public synchronized boolean removeAction(Action oldAction) {
        if (super.removeAction(oldAction)) {
            actionStatusMap.remove(oldAction);
            actionExceptionMap.remove(oldAction);
            return true;
        }
        return false;
    }
    
    /**
     * Returns the status of an action held by the controller.  Returns
     * null for actions not held by the controller.
     */
    public synchronized ActionStatus getActionStatus(Action action) {
        return (ActionStatus)actionStatusMap.get(action);
    }
    
    /**
     * Returns the exception thrown by an action held by the controller.  
     * Returns null if the given action did not throw an exception or is not 
     * held by the controller.
     */
    public synchronized Exception getException(Action action) {
        return (Exception)actionStatusMap.get(action);
    }
    
    /**
     * Causes uncompleted actions added to this group to be run in sequence.
     * Should not be executed directly, but instead it is executed by a thread
     * made upon invoking the start method.
     */
    public void run() {
        Thread.currentThread().setName("Controller thread");
        fireEvent(new ControllerEvent(this, ControllerEvent.START));
        while(numActions > 0) {
            synchronized(this) {
                currentAction = pendingActions[0];
                //removeAction will remove the action from the status map
                removeAction(currentAction);
                //put it back
                actionStatusMap.put(currentAction,CURRENT);
            }
            fireEvent(new ControllerEvent(this, ControllerEvent.START_ACTION, currentAction));
            if(currentAction instanceof Activity) {
                waitObject.currentActionDone = false;
                waitObject.actionException = null;
                //define a thread to run the activity
                Thread localRunner = new Thread(new Runnable() {
                    public void run() {
                        Thread.currentThread().setName(currentAction.getLabel()+" thread");
                        try {
                            currentAction.actionPerformed();
                        }
                        catch (RuntimeException e) {
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
                            //need to check again that current action isn't done, to ensure that current-action thread didn't finish between while and here
                            while(urgentAction == null && !waitObject.currentActionDone) {
                                waitObject.wait();
                            }
                        }
                    } catch(InterruptedException e) {}

                    if (waitObject.actionException == null) {
                        synchronized(this) {
                            doUrgentAction();
                            if(!wasPaused) unPause();
                        }
                    }
                }
            } else {//currentAction is not an Activity; run on group's thread
                currentAction.actionPerformed();
            }
            
            synchronized(this) {
                //update the action's status
                actionStatusMap.remove(currentAction);
                if (waitObject.actionException != null) {
                    actionStatusMap.put(currentAction,FAILED);
                    actionExceptionMap.put(currentAction,waitObject.actionException);
                }
                else if (haltRequested) {
                    actionStatusMap.put(currentAction,STOPPED);
                }
                else {
                    actionStatusMap.put(currentAction,COMPLETED);
                }
                completedActions = (Action[])Arrays.addObject(completedActions, currentAction);
                fireEvent(new ControllerEvent(this, ControllerEvent.END_ACTION, currentAction));

                if(waitObject.actionException == null && repeatCurrentAction) {
                    if(currentAction instanceof Activity) {
                        addNextAction(((Activity)currentAction).makeCopy());
                    } else {
                        addNextAction(currentAction);
                    }
                }
            
                currentAction = null;
            }

            if (waitObject.actionException != null) {
                fireEvent(new ControllerEvent(this, ControllerEvent.NO_MORE_ACTIONS));
                throw new RuntimeException("action failed", waitObject.actionException);
            }
            
            doUrgentAction();
            
            if(pauseAfterEachAction) {
                pauseRequested = true;
            }
            
            if(haltRequested) {
                fireEvent(new ControllerEvent(this, ControllerEvent.HALTED));
                haltRequested = false;
            }
            
            if(!doContinue()) {
                break;
            }

        }//end while(numActions > 0)
        doUrgentAction();//could come straight here if doActionNow before adding other actions
        synchronized(this) {
            notifyAll();//notify any threads requesting halt and waiting for execution to complete
        }

        fireEvent(new ControllerEvent(this, ControllerEvent.NO_MORE_ACTIONS));
    }

    /**
     * Pauses current activity, executes given action, then resumes current
     * activity. If current action is already paused, it is not resumed by this
     * method.
     * 
     * @param action
     *            Action to be performed right away; cannot be an Activity.
     */
    public synchronized void doActionNow(final Action action) {
        if(action instanceof Activity) throw new IllegalArgumentException("can't send Activity here; sorry");
//        System.out.println("in doActionNow, urgentAction="+urgentAction);
        if(urgentAction != null) return;
        if(isActive()) {
            wasPaused = isPaused();
//            System.out.println("in doActionNow "+wasPaused);
            pause();
//            System.out.println("in doActionNow paused it");
            synchronized(waitObject) {
                urgentAction = action;
                waitObject.notifyAll();//release group thread from wait
            }
        } else {
            Thread localRunner = new Thread(new Runnable() {
                public void run() {
                    synchronized(Controller.this) {
                        if(Controller.this.isActive()) {//check this just in case something calls start between time when gui thread exits method and when localRunner.run() is started
                            doActionNow(action);
                        } else {
                            action.actionPerformed();
                        }
                    }
                }
            });
            localRunner.start();
        }
    }//end of run
    
    private synchronized void doUrgentAction() {
//        System.out.println("doing UrgentAction "+urgentAction);
        if(urgentAction == null) return;
        fireEvent(new ControllerEvent(this, ControllerEvent.START_URGENT_ACTION, urgentAction));
        urgentAction.actionPerformed();
        completedActions = (Action[])Arrays.addObject(completedActions, urgentAction);
//        System.out.println("finished UrgentAction "+urgentAction);
        fireEvent(new ControllerEvent(this, ControllerEvent.END_URGENT_ACTION, urgentAction));
        urgentAction = null;
    }
    
    /**
     * Adds an action to the beginning of the pending actions list, so that it will be performed next.
     * This is used by the repeatCurrentAction facility.
     */
    private synchronized void addNextAction(Action nextAction) {
        Action[] newActions = new Action[pendingActions.length+1];
        newActions[0] = nextAction;
        for(int i=0; i<pendingActions.length; i++) {
            newActions[i+1] = pendingActions[i];
        }
        pendingActions = newActions;
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
    
    /**
     * Request that the current activity terminate and wait till it does.
     * Controller is left in a paused state after termination of activity.
     * If current action is not an Activity, this has the same effect as pause.  
     * Has no effect on Controller itself, or any pending activities, 
     * which remain in the queue.
     */
    public synchronized void halt() {
        if(!isActive()) return;
        pauseRequested = true;
        haltRequested = true;
        if(currentAction instanceof Activity) {
            ((Activity)currentAction).halt();
        }
    }

    /**
     * Adds the given object to the list of listeners to this Controller.
     * Nothing is done if listener is null.
     */
    public synchronized void addListener(ControllerListener listener) {
        if(listener == null) return;
        eventManager.addListener(listener);
    }

    /**
     * Removes the given object from the list of listeners to this Controller.
     */
    public synchronized void removeListener(ControllerListener listener) {
        if(listener == null) return;
        eventManager.removeListener(listener);
    }

    /**
     * Notifies all registered listeners of the given event.
     */
    protected synchronized void fireEvent(ControllerEvent event) {
        eventManager.fireEvent(event);
    }
    
    public String toString() {
        return "Controller";
    }
    
    private final ControllerEventManager eventManager;

    private boolean wasPaused = false;
    private Action urgentAction;
    
    private boolean repeatCurrentAction = false;

    protected final WaitObject waitObject;
    private static class WaitObject implements java.io.Serializable {
        boolean currentActionDone;
        Exception actionException;
    };

    private final HashMap actionStatusMap, actionExceptionMap;
    
    /**
     * Enumerated type describing the status of an action.
     */
    public static class ActionStatus extends EnumeratedType {
        protected ActionStatus(String label) {
            super(label);
        }       
        public EnumeratedType[] choices() {return CHOICES;}
        /**
         * Required to guarantee singleton when deserializing.
         * @return the singleton INSTANCE
         */
        private Object readResolve() {
            for (int i=0; i<CHOICES.length; i++) {
                if (this.toString().equals(CHOICES[i].toString())) {
                    return CHOICES[i];
                }
            }
            throw new RuntimeException("unknown action status: "+this);
        }
    }
    protected static final ActionStatus[] CHOICES = 
        new ActionStatus[] {
            new ActionStatus("Pending"),
            new ActionStatus("Current"), 
            new ActionStatus("Completed"),
            new ActionStatus("Stopped"),
            new ActionStatus("Failed")};
    public static final ActionStatus PENDING = CHOICES[0];
    public static final ActionStatus CURRENT = CHOICES[1];
    public static final ActionStatus COMPLETED = CHOICES[2];
    public static final ActionStatus STOPPED = CHOICES[3];
    public static final ActionStatus FAILED = CHOICES[4];

}
