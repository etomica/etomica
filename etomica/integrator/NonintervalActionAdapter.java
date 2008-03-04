package etomica.integrator;

import etomica.api.IAction;
import etomica.api.IIntegratorNonintervalListener;
import etomica.integrator.IntegratorNonintervalEvent.NonintervalEventType;

/**
 * Adapter that causes an action to be performed as the result of an integrator
 * non-interval event. Types of events that trigger action can be set using the
 * setEventTypes method; default is to respond to all events.
 */
public class NonintervalActionAdapter implements IIntegratorNonintervalListener, java.io.Serializable {

    /**
     * Creates adapter with integrator to be set later.
     */
    public NonintervalActionAdapter(IAction action) {
        this.action = action;
        setActive(true);
        setEventTypes((NonintervalEventType[])NonintervalEventType.choices());
    }

    /**
     * Method of Integrator.IntervalListener interface. After receiving an
     * appropriate event updateInterval times, the action will be performed. If
     * event is not of the type indicated by setEventTypes, no action is
     * performed and the updateInterval counter is not incremented; this is also
     * the case if the active flag is false.
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt) {
        if (active) {
            for (int i=0; i<eventTypes.length; i++) {
                if (eventTypes[i] == evt.type()) {
                    action.actionPerformed();
                    return;
                }
            }
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
    public IAction getAction() {
        return action;
    }

    /**
     * Sets the types of integrator event that triggers action.
     * 
     * @param types
     *            array of types that will trigger the action
     */
    public void setEventTypes(NonintervalEventType[] types) {
        eventTypes = (NonintervalEventType[])types.clone();
    }

    private static final long serialVersionUID = 1L;
    private final IAction action;
    private NonintervalEventType[] eventTypes;
    private boolean active;//set true in constructor
}
