package etomica;

/**
 * Interface for an object that can be registered as a listener to another simulation element.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface MCMoveEventListener extends SimulationEventListener {

    public void mcMoveAction(MCMoveEvent event);

}
