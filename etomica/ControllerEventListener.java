package etomica;

/**
 * Interface for an object that can be registered as a listener to another simulation element.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface ControllerEventListener extends SimulationEventListener {

    public void controllerAction(ControllerEvent event);

}
