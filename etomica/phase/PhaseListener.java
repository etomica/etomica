package etomica.phase;


/**
 * Interface for an object that can be registered as a listener to another simulation element.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface PhaseListener {

    public void actionPerformed(PhaseEvent event);

}
