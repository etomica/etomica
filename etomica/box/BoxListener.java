package etomica.box;


/**
 * Interface for an object that can be registered as a listener to another simulation element.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface BoxListener {

    public void actionPerformed(BoxEvent event);

}
