package etomica.integrator;

import etomica.action.Action;

/**
 * Adapter that causes an action to be performed as the result of an integrator
 * non-interval event. Types of events that trigger action can be set using the
 * setEventTypes method; default is to respond to all events.
 */

/*
 * History Created on Feb 19, 2005 by kofke
 */
public class NonintervalActionAdapter implements IntegratorNonintervalListener, java.io.Serializable {

    /**
     * Creates adapter such that given action is performed in response to
     * interval events from the given integrator.
     */
    public NonintervalActionAdapter(Action action, Integrator integrator) {
        this(action);
        integrator.addListener(this);
    }

    /**
     * Creates adapter with integrator to be set later.
     */
    public NonintervalActionAdapter(Action action) {
        this.action = action;
        setActive(true);
        setEventTypes(IntegratorNonintervalEvent.Type.choices);
    }

    /**
     * Method of Integrator.IntervalListener interface. After receiving an
     * appropriate event updateInterval times, the action will be performed. If
     * event is not of the type indicated by setEventTypes, no action is
     * performed and the updateInterval counter is not incremented; this is also
     * the case if the active flag is false.
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (active && ((evt.type().mask & eventMask) != 0)) {
            action.actionPerformed();
        }
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
     * @return Returns the action.
     */
    public Action getAction() {
        return action;
    }

    /**
     * Sets the types of integrator event that triggers action.
     * 
     * @param types
     *            array of types that will trigger the action
     */
    public void setEventTypes(IntegratorIntervalEvent.Type[] types) {
        eventMask = 0;
        for (int i = 0; i < types.length; i++) {
            eventMask |= types[i].mask;
        }
    }

    private final Action action;
    private int eventMask;
    private boolean active;//set true in constructor

}
