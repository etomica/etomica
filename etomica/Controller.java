package etomica;

import java.util.LinkedList;

/**
 * Organizer of actions of the integrators.  Controller sets the protocol for
 * conduct of the simulation.  For example, the simulation might start and stop 
 * with the push of a button; it might run for a fixed number of cycles and exit;
 * it might run several simulations over a series of states.<p>
 * The Controller runs on its own thread, and it spawns Integrator processes
 * each on its own thread.
 */
public class Controller extends ActivityGroupSeries implements Runnable, java.io.Serializable, EtomicaElement {
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Simple controller that enables start/stop of simulation using a button");
        return info;
    }
 
    /**
     * Method to start this controller's thread, causing execution of the run() method.
     * No action is performed if a thread is already running this controller (i.e., if
     * start was called before, and actions remain to be completed).
     */
    public void start() {
        if(isActive()) return;
        pauseRequested = false;
        runner = new Thread(this);
        runner.start();
    }

    public synchronized void addListener(ControllerListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(ControllerListener listener) {
        eventManager.removeListener(listener);
    }

    protected void fireEvent(ControllerEvent event) {
        eventManager.fireEvent(event);
    }    
    
    private SimulationEventManager eventManager = new SimulationEventManager();

}//end of Controller


