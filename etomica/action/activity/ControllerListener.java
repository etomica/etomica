package etomica.action.activity;


/**
 * Interface for an object that can be registered as a listener to another simulation element.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface ControllerListener {

    public void actionPerformed(ControllerEvent event);

}
