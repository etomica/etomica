package etomica.integrator.mcmove;

import etomica.SimulationListener;

/**
 * Interface for an object that listens for events fired by a MCMove object.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface MCMoveListener extends SimulationListener {

    public void actionPerformed(MCMoveEvent event);

}
