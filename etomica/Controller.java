package etomica;

//Java2 imports
import java.util.LinkedList;
import java.util.Iterator;

//import etomica.utility.java2.LinkedList;
//import etomica.utility.java2.Iterator;

/**
 * Organizer of activities of the integrators.  Controller sets the protocol for
 * conduct of the simulation.  For example, the simulation might start and stop 
 * with the push of a button; it might run for a fixed number of cycles and exit;
 * it might run several simulations over a series of states.<p>
 * The Controller runs on its own thread, and it spawns Integrator processes
 * each on its own thread.
 */
public class Controller implements Runnable, java.io.Serializable, EtomicaElement {

  /**
   * List of integrators managed by the controller
   */
    protected final LinkedList integrators = new LinkedList();
  /**
   * Thread used to run the controller
   */
    protected transient Thread runner;
    
    private SimulationEventManager eventManager = new SimulationEventManager();

    /**
     * Creates a new controller.
     *
     */
    public Controller() {
        maxSteps = Integer.MAX_VALUE;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Simple controller that enables start/stop of simulation using a button");
        return info;
    }

   /**
    * Resets the controller in some manner.
    * Performs no action for this controller, but may do something in subclasses
    */
    public void reset() {}
 
   /**
    * Adds an integrator to the list of integrators managed by this controller.
    * Starts the integrator immediately if autoStart is true
    */
    public void add(Integrator i) {
        integrators.add(i);
        i.setController(this);
        i.setMaxSteps(maxSteps);
        if(autoStart) i.start();
    }
    
    /**
     * @return a list of the integrators managed by this controller
     */
    public LinkedList integrators() {return integrators;}
    
    /**
     * Accessor method for the autoStart flag.
     * autoStart flags for whether the integrators start as soon as added to the controller.
     * Default is false
     */
    public void setAutoStart(boolean b) {autoStart = b;}
    /**
     * Accessor method for the autoStart flag.
     * autoStart flags for whether the integrators start as soon as added to the controller.
     * Default is false
     */
    public boolean getAutoStart() {return autoStart;}
    
    /**
     * Accessor for maxSteps, the total number of steps performed by each integrator.
     */
    public int getMaxSteps() {return maxSteps;}
    /**
     * Accessor for maxSteps, the total number of steps performed by each integrator.
     */
    public void setMaxSteps(int m) {
        maxSteps = m;
        for(Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            ((Integrator)iter.next()).setMaxSteps(m);
        }
    }

    /**
     * Method to start this controller's thread, causing execution of the run() method
     */
    public void start() {
        if(runner != null) return;
    	active = true;
        paused = false;
        runner = new Thread(this);
        runner.start();
    }
                    
    /**
     * Activity performed by this controller once its start method is called.
     * This controller simply initiates each integrator on its own thread, which
     * is initiated by calling its start method.  Controller subclasses would override
     * this method to perform some other function.<p>
     * An integrator can be run on the controller's thread by instead calling
     * the run() method of the integrator.
     */
    public void run() {
        for(Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            if(!integrator.isActive()) integrator.start();
        }
        runner = null;
    }

    public void runInDevelopment() {
    	LinkedList list;
//    	while(true) {
//    		//get next SimulationActivity
//    		//activity
//    		//activity.join
//    	}
        for(Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            if(!integrator.isActive()) integrator.start();
        }
        runner = null;
    }

    /**
     * Sends a halt (stop) signal to all integrators.
     */
    public void halt() {
    	active = false;
        for(Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            if(integrator.isActive()) integrator.halt();
        }
        paused = true;
    }
    
    /**
     * Sends a pause signal to all active integrators.
     */
    public void pause() {
        if(paused) return;
        paused = true;
        for(Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            if(integrator.isActive()) integrator.pause();
        }
    }

    /**
     * Sends a resume signal to all active integrators.
     */
    public void unPause() {
        if(!paused) return;
        paused = false;
        for(Iterator iter=integrators.iterator(); iter.hasNext(); ) {
            Integrator integrator = (Integrator)iter.next();
            if(integrator.isActive()) integrator.unPause();
        }
    }
    
    /**
     * Returns true if any integrator governed by this controller is active.
     * Returns false if none are active.
     */
    public boolean isActive() {
    	return active;
//        for(Iterator iter=integrators.iterator(); iter.hasNext(); ) {
//            if(((Integrator)iter.next()).isActive()) return true;
//        }
//        return false;
    }
    
    /**
     * Returns true if controller is not active or if pause() was last called
     * without an intervening call to unPause().
     */
    public boolean isPaused() {return paused;}
        
    
    public synchronized void addListener(ControllerListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(ControllerListener listener) {
        eventManager.removeListener(listener);
    }

    protected void fireEvent(ControllerEvent event) {
        eventManager.fireEvent(event);
    }    
    
    private boolean initialized = false;
    private boolean autoStart = false;
    private int maxSteps;
    private boolean paused = true;
    private boolean active = false;

}//end of Controller


