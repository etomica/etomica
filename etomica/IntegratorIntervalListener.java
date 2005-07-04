package etomica;

/**
 * Interface for a class that listens in particular to interval events from an
 * integrator.  Interval events are fired regularly by the integrator to announce
 * its progress with the simulation.  In contrast, non-interval events are fired
 * at special points in the simulation, announcing for example the beginning and
 * end of the simulation.
 */

public interface IntegratorIntervalListener extends IntegratorListener {

    /**
     * Action performed by the listener when integrator fires its interval
     * event.
     */
    public void intervalAction(IntegratorIntervalEvent evt);

    /**
     * Priority assigned to the listener. A small value will cause the
     * listener's action to be performed earlier, before another listener having
     * a larger priority value (e.g., priority 100 action is performed before
     * one having a priority of 200). Listeners that cause periodic boundaries
     * to be applied are given priorities in the range 100-199. Ordering is
     * performed only when a listener is added to the integrator, or when the
     * sortListeners method of integrator is called.
     */
    public int getPriority();
}