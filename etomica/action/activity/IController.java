package etomica.action.activity;

import etomica.action.Action;
import etomica.action.activity.Controller.ActionStatus;

public interface IController extends ActivityGroup {

    /**
     * Marks all actions as pending.  Previously completed actions are put 
     * first in the new list of pending actions.  The Controller must be 
     * inactive or halted before calling this method.
     */
    public void reset();

    /**
     * Returns the status of an action held by the controller.  Returns
     * null for actions not held by the controller.
     */
    public ActionStatus getActionStatus(Action action);

    /**
     * Returns the exception thrown by an action held by the controller.  
     * Returns null if the given action did not throw an exception or is not 
     * held by the controller.
     */
    public Exception getException(Action action);

    /**
     * Pauses current activity, executes given action, then resumes current
     * activity. If current action is already paused, it is not resumed by this
     * method.
     * 
     * @param action
     *            Action to be performed right away; cannot be an Activity.
     */
    public void doActionNow(final Action action);

    /**
     * Returns the event manager used by the controller to notify listeners of
     * individual actions events as well as events related to the controller
     * itself.
     * @see ControllerEvent
     */
    public ControllerEventManager getEventManager();

}