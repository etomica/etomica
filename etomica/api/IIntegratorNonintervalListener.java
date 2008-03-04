package etomica.api;

import etomica.integrator.IntegratorNonintervalEvent;


/**
 * Interface for a class that listens for integrator events other than interval events.
 * Non-interval events are thrown special points in the simulation, such as when it
 * begins and ends.
 */
public interface IIntegratorNonintervalListener {

    /**
     * Action performed by the listener when integrator fires its event.
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt);

}
