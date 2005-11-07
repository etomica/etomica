package etomica.integrator.mcmove;


/**
 * Interface for an object that listens for events fired by a MCMove object.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface MCMoveListener {

    public void actionPerformed(MCMoveEvent event);

}
