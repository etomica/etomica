package etomica;

/**
 * Interface for an object that can be registered as a listener to another simulation element.
 * Can be managed by a SimulationEventManager.
 *
 * @author David Kofke
 */
public interface SimulationEventListener extends java.util.EventListener {

//    public void simulationAction(SimulationEvent event);

    /**
     * Class used to construct a two-way linked list of listeners.
     */
    public static class Linker implements java.io.Serializable {
        public SimulationEventListener listener;
        public Linker next, previous;
        public Linker(SimulationEventListener listener, Linker next, Linker previous) {
            this.listener = listener;
            this.next = next;
            this.previous = previous;
        }
    }
}
