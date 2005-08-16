package etomica;

import etomica.action.activity.ActivityGroupSeries;
import etomica.utility.Arrays;

/**
 * Organizer of actions of the integrators.  Controller sets the protocol for
 * conduct of the simulation.  For example, the simulation might start and stop 
 * with the push of a button; it might run for a fixed number of cycles and exit;
 * it might run several simulations over a series of states.<p>
 * The Controller runs on its own thread, and it spawns Integrator processes
 * each on its own thread.
 */
public class Controller extends ActivityGroupSeries implements java.io.Serializable, EtomicaElement {
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Simple controller that enables start/stop of simulation using a button");
        return info;
    }
 
    /**
     * Method to start this controller's thread, causing execution of the run() method.
     * No action is performed if a thread is already running this controller (i.e., if
     * start was called before, and actions remain to be completed).
     */
    public synchronized void start() {
        if(isActive()) return;
        runner = new Thread(new Runnable() {
            public void run() {actionPerformed();}
        });
        runner.start();
    }
    
    /**
     * Causes uncompleted actions added to this group to be run in sequence.  Should not be
     * executed directly, but instead it is executed
     * by a thread made upon invoking the start method.
     */
    public void run() {
        Thread.currentThread().setName("Controller thread");
        fireEvent(new ControllerEvent(this, ControllerEvent.START));
        while(numActions > 0) {
            synchronized(this) {
                currentAction = actions[0];
                removeAction(currentAction);
            }
            fireEvent(new ControllerEvent(this, ControllerEvent.START_ACTION, currentAction));
            if(currentAction instanceof Activity) {
                waitObject.currentActionDone = false;
                waitObject.exceptionThrown = false;
                //define a thread to run the activity
                Thread localRunner = new Thread(new Runnable() {
                    public void run() {
                        Thread.currentThread().setName(currentAction.getLabel()+" thread");
                        RuntimeException exception = null;
                        try {
                            currentAction.actionPerformed();
                        }
                        catch (RuntimeException e) {
                            exception = e;
                            synchronized(waitObject) {
                                waitObject.exceptionThrown = true;
                            }
                        }
                        synchronized(waitObject) {
                            waitObject.currentActionDone = true;
                            waitObject.notifyAll();
                        }
                        if (exception != null) {
                            throw exception;
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
//                            System.out.println("done waiting for waitObject");
                        }
                    } catch(InterruptedException e) {}

                    if (waitObject.exceptionThrown) {
                        throw new RuntimeException("action failed");
                    }
                    
                    synchronized(this) {
//                        System.out.println("trying to do urgentAction");
                        doUrgentAction();
//                        System.out.println("in run "+wasPaused);
                        if(!wasPaused) unPause();
                    }
                    
                }
            } else {//currentAction is not an Activity; run on group's thread
                currentAction.actionPerformed();
            }
            
            //TODO mark this as whether completed normally
            synchronized(this) {
                completedActions = (Action[])Arrays.addObject(completedActions, currentAction);
                fireEvent(new ControllerEvent(this, ControllerEvent.END_ACTION, currentAction));

                if(repeatCurrentAction) {
                    if(currentAction instanceof Activity) {
                        addNextAction(((Activity)currentAction).makeCopy());
                    } else {
                        addNextAction(currentAction);
                    }
                }
            
                currentAction = null;
            }

            doUrgentAction();
            
            if(pauseAfterEachAction) {
                pauseRequested = true;
            }
            
            if(!doContinue()) {
                break;
            }

        }//end while(numActions > 0)
        doUrgentAction();//could come straight here if doActionNow before adding other actions
        synchronized(this) {
            notifyAll();//notify any threads requesting halt and waiting for execution to complete
        }
        if(haltRequested) {
            fireEvent(new ControllerEvent(this, ControllerEvent.HALTED));
        } else {
            fireEvent(new ControllerEvent(this, ControllerEvent.NO_MORE_ACTIONS));
        }
    }

    /**
     * Pauses current activity, executes given action, then resumes current activity.
     * If current action is already paused, it is not resumed by this method.
     * @param action Action to be performed right away; cannot be an Activity.
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
    }
    
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
        actions = (Action[])Arrays.resizeArray(actions, actions.length+1);
        if(actions.length > 1) {
            System.arraycopy(actions, 0, actions, 1, actions.length-1);
        }
        actions[0] = nextAction;
    }
    
    /**
     * Flag indicating whether the current action should be repeated indefinitely before
     * moving on to the next action.  Setting this to true also sets pauseAfterEachAction 
     * to true (setting false then has no effect on pauseAfterEachAction).  In repeating
     * an action, when an action is completed it is placed again at the head of the pending 
     * list (for an Action the same instance is use; for an Activity a copy is made).
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
    
    public synchronized boolean isActive() {
    	    return (runner != null) && runner.isAlive();
    }
    
    /**
     * Request that the current activity terminate and wait till it does.  
     * Has no effect on Controller itself, or any pending activities, 
     * which remain in the queue.
     */
    public synchronized void halt() {
        if(!isActive()) return;
        pauseRequested = true;
        if(currentAction instanceof Activity) {
            ((Activity)currentAction).halt();
        }
    }


    public synchronized void addListener(ControllerListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(ControllerListener listener) {
        eventManager.removeListener(listener);
    }

    protected synchronized void fireEvent(ControllerEvent event) {
        eventManager.fireEvent(event);
    }    
    
    private SimulationEventManager eventManager = new SimulationEventManager();

    private Thread runner;
    
    private boolean wasPaused = false;
    private Action urgentAction;
    
    private boolean repeatCurrentAction;

    private final WaitObject waitObject = new WaitObject();
    private static class WaitObject implements java.io.Serializable {
        boolean currentActionDone;
        boolean exceptionThrown;
    };

}//end of Controller


