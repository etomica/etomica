package etomica.lattice;
import etomica.SimulationListener;

/**
 * Interface for an object that can be registered as a listener to a lattice.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface LatticeListener extends SimulationListener {

    public void actionPerformed(LatticeEvent event);

}
