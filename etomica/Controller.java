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
        while(numActions > 0) {
            synchronized(this) {
                currentAction = actions[0];
                removeAction(currentAction);
            }

            if(currentAction instanceof Activity) {
                waitObject.currentActionDone = false;
                waitObject.exceptionThrown = false;
                Thread localRunner = new Thread(new Runnable() {
                    public void run() {
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
                localRunner.start();
                while(!waitObject.currentActionDone) {
                    try {
                        synchronized(waitObject) {
                            //put group thread in wait state while current action runs
                            //need to check again that current action isn't done, to ensure that current-action thread didn't finish between while and here
                            if(urgentAction == null && !waitObject.currentActionDone) {
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
                currentAction = null;
            }
            doUrgentAction();
            if(pauseRequested || pauseAfterEachAction) doWait();
            if(haltRequested) break;
        }
        doUrgentAction();//could come straight here if doActionNow before adding other actions
        synchronized(this) {
            notifyAll();//notify any threads requesting halt and waiting for execution to complete
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
        urgentAction.actionPerformed();
        completedActions = (Action[])Arrays.addObject(completedActions, urgentAction);
//        System.out.println("finished UrgentAction "+urgentAction);
        urgentAction = null;
    }
    
    public synchronized boolean isActive() {
    	return (runner != null) && runner.isAlive();
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
    
    private SimulationEventManager eventManager = new SimulationEventManager();

    private Thread runner;
    
    private boolean wasPaused = false;
    private Action urgentAction;

    private final WaitObject waitObject = new WaitObject();
    private static class WaitObject {
        boolean currentActionDone;
        boolean exceptionThrown;
    };

}//end of Controller


