package etomica.integrator;

import etomica.Action;
import etomica.Integrator;
import etomica.IntegratorIntervalEvent;
import etomica.IntegratorIntervalListener;

/**
 * Adapter that causes an action to be performed as the result of an integrator
 * intervalEvent.
 */

/*
 * History Created on Feb 19, 2005 by kofke
 */
public class IntervalActionAdapter implements IntegratorIntervalListener {

    /**
     * Creates adapter such that given action is performed in response to
     * interval events from the given integrator.
     */
    public IntervalActionAdapter(Action action, Integrator integrator) {
        this(action);
        integrator.addIntervalListener(this);
    }

    /**
     * Creates adapter with integrator to be set later.
     */
    public IntervalActionAdapter(Action action) {
        this.action = action;
        setActive(true);
        setActionInterval(1);
        setPriority(400);
    }

    /**
     * Method of Integrator.IntervalListener interface. After receiving an
     * appropriate event updateInterval times, the action will be performed. If
     * event is not of the type indicated by setEventTypes, no action is
     * performed and the updateInterval counter is not incremented; this is also
     * the case if the active flag is false.
     */
    public void intervalAction(IntegratorIntervalEvent evt) {
        if (active) {
            if (--iieCount == 0) {
                iieCount = actionInterval;
                action.actionPerformed();
            }
        }
    }

    /**
     * Returns the priority, which determines the order in which listeners are
     * informed when integrator fires an interval event. Default is 400.
     * 
     * @see etomica.IntegratorIntervalListener#getPriority()
     */
    public int getPriority() {
        return priority;
    }

    /**
     * Sets the priority, which determines the order in which listeners are
     * informed when integrator fires an interval event. Default is 400.
     */
    public void setPriority(int priority) {
        this.priority = priority;
    }

    /**
     * @return Flag indicating whether instance responds to interval events
     */
    public boolean isActive() {
        return active;
    }

    /**
     * Sets flag indicating whether instance responds to interval events.
     * Default is true.
     */
    public void setActive(boolean active) {
        this.active = active;
    }

    /**
     * Returns the number of interval events between each action performed.
     */
    public int getActionInterval() {
        return actionInterval;
    }

    /**
     * Sets the number of interval events between each action performed. Default
     * is 1.
     */
    public void setActionInterval(int actionInterval) {
        this.actionInterval = actionInterval;
        iieCount = actionInterval;
    }

    /**
     * @return Returns the action.
     */
    public Action getAction() {
        return action;
    }

    private final Action action;
    private int iieCount;
    private int priority;
    private boolean active;//set true in constructor
    private int actionInterval;

}