package etomica;


/**
 * Interface for a class that listens for integrator events, typically
 * those other than interval event.
 *
 */

/*
 * History
 * Created on Apr 26, 2005 by kofke
 */
public interface IntegratorNonintervalListener extends IntegratorListener {

    /**
     * Action performed by the listener when integrator fires its event.
     */
    public void nonintervalAction(IntegratorNonintervalEvent evt);

}
