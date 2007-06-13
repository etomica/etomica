package etomica.integrator;


/**
 * Interface for a class that listens for integrator events other than interval events.
 * Non-interval events are thrown special points in the simulation, such as when it
 * begins and ends.
 */
public interface IntegratorNonintervalListener {

    /**
     * Action performed by the listener when integrator fires its event.
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt);

}
