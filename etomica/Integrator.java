package etomica;
import java.util.*;
import etomica.units.*;

/**
 * Integrator is used to define the algorithm used to move the atoms around and
 * generate new configurations in one or more phases.  All integrator methods, 
 * such as molecular dynamics or Monte Carlo are implemented via subclasses of
 * this Integrator class.  The Integrator's activities are managed via the 
 * actions of the governing Controller.
 *
 * @author David Kofke
 */
public abstract class Integrator extends SimulationElement implements Runnable, java.io.Serializable {

  public static String VERSION = "Integrator:01.12.04";
  
  transient public Thread runner;
  private boolean haltRequested = false;
  private boolean resetRequested = false;
  protected int maxSteps = Integer.MAX_VALUE;
  protected int stepCount = 0;
  //looking to eliminate isPaused field and used just pauseRequested
//  private boolean isPaused = true;
  protected boolean pauseRequested = false;
  protected final PotentialMaster potential;
 
  protected Phase firstPhase;
  protected Phase[] phase;
  int phaseCount = 0;
  int phaseCountMax = 1;
  protected int sleepPeriod = 10;
  
  //should use a different collection structure
  private Vector intervalListenersBeforePbc = new Vector();
  private Vector intervalListenersImposePbc = new Vector();
  private Vector intervalListenersAfterPbc = new Vector();
  int interval = 10;  // number of steps between IntervalEvent firing
  int integrationCount = 0;
  boolean doSleep = true;
  protected double temperature = Default.TEMPERATURE;
  protected boolean isothermal = false;
  private boolean initialized = false;
  private boolean running = false;
  
  IntervalEvent intervalEvent = new IntervalEvent(this, IntervalEvent.INTERVAL);
  public Controller parentController;

  public Integrator(Simulation sim) {
    super(sim, Integrator.class);
    phase = new Phase[phaseCountMax];
    potential = sim.hamiltonian.potential;
  }
  
    /** 
     * Performs the elementary integration step, such as a molecular dynamics
     * time step, or a Monte Carlo trial.
     */
    public abstract void doStep();
    
    /**
     * Defines the actions taken by the integrator to reset itself, such
     * as required if a perturbation is applied to the simulated phase by (e.g., another
     * thread adds or deletes a molecule).  Also invoked when the integrator
     * is started or initialized.  This method is called by the integrator's thread
     * when running if the reset() method was previously called on another thread.
     * This (doReset) method should not be called by another class.
     */
    protected abstract void doReset(); //protected because runs on integrator's thread
    
    /**
     * Returns a new instance of an agent of this integrator for placement in
     * the given atom in the ia (IntegratorAgent) field.
     */
    public abstract Agent makeAgent(Atom a);
  
  /**
   * Initializes the integrator, performing the following steps:  (1) deploys
   * agents in all atoms; (2) call doReset method; (3) fires an event indicating
   * to registered listeners indicating that initialization has been performed
   * (i.e. fires IntervalEvent of type field set to INITIALIZE).
   */
  public void initialize() {
    deployAgents();
    initialized = true;
    doReset();
    fireIntervalEvent(new IntervalEvent(this, IntervalEvent.INITIALIZE));
  }
  /**
   * Returns true if the initialize method has been called, and the 
   * run method has not completed.
   */
  public boolean isInitialized() {return initialized;}
  
  public void setController(Controller c) {parentController = c;}
  
  //how do agents get placed in atoms made during the simulation?
  protected void deployAgents() {  //puts an Agent of this integrator in each atom of all phases
    for(int i=0; i<phaseCount; i++) {
        Phase p = phase[i];
        for(Atom a=p.firstAtom(); a!=null; a=a.nextAtom()) {
            a.setIntegratorAgent(makeAgent(a));
        }
    }
  }
  
  public final void setTemperature(double t) {temperature = t;}
  public final double getTemperature() {return temperature;}
  public final double temperature() {return temperature;}
  public Dimension getTemperatureDimension() {return Dimension.TEMPERATURE;}

  //Other introspected properties
  public void setIsothermal(boolean b) {isothermal = b;}
  public boolean isIsothermal() {return isothermal;}
    
  public final int getInterval() {return interval;}
  public final void setInterval(int interval) {this.interval = interval;}
  
  public final int getSleepPeriod() {return sleepPeriod;}
  public final void setSleepPeriod(int s) {sleepPeriod = s;}

  public void setDoSleep(boolean b) {doSleep = b;}
  public boolean isDoSleep() {return doSleep;}

  /**
   * @return true if integrator can perform integration of another phase, 
   *         false if the integrator has all the phases it was built to handle
   */
  public boolean wantsPhase() {return phaseCount < phaseCountMax;}
  
  /**
   * Calls the setPhase method with the given phase.
   * @deprecated  instead use setIntegrator method in phase.
   */
  public void registerPhase(Phase p) {
    setPhase(p);
  }
  
  /**
   * @deprecated  instead use setIntegrator method in phase.
   */
  public void setPhase(Phase p) {
    addPhase(p);
  }
  
  public Phase getPhase(int i) {return phase[i];}
  
  /**
   * Performs activities needed to set up integrator to work on given phase.
   * This method should not be called directly; instead it is invoked by the phase in its setIntegrator method.
   * @return true if the phase was successfully added to the integrator; false otherwise
   */
   //perhaps should throw an exception rather than returning a boolean "false"
  public boolean addPhase(Phase p) {
    for(int i=0; i<phaseCount; i++) {if(phase[i]==p) return false;}  //check that phase is not already registered
    if(!this.wantsPhase()) {return false;}  //if another phase not wanted, return false
    phase[phaseCount] = p;
    phaseCount++;
    firstPhase = phase[0];
    makeIterators(p.iteratorFactory());
	p.iteratorFactoryMonitor.addObserver(iteratorFactoryObserver());
    return true;
  }
  
  /**
   * Performs activities needed to disconnect integrator from given phase.
   * This method should not be called directly; instead it is invoked by the phase in its setIntegrator method
   */
  public void removePhase(Phase p) {
    for(int i=0; i<phaseCount; i++) {
        if(phase[i]==p) {//phase found; remove it
            phase[i] = null;
            phaseCount--;
            if(phaseCount > 0) phase[i] = phase[phaseCount];
            firstPhase = phase[0];
            break;
        }
    }
  }
  
  /**
   * Attempts to set the given integrator as the integrator for all phases
   * that were added to this phase.
   */
  public void transferPhasesTo(Integrator anotherIntegrator) {
    for(int i=0; i<phaseCount; i++) {
        phase[i].setIntegrator(anotherIntegrator);
    }
  }

	/**
	 * Method performs no action, but can be overridden in subclasses to handle setting or change of iteratorFactory in phase
	 */
	 //this should be made abstract
	protected void makeIterators(IteratorFactory i) {}

	/**
	 * Returns an observer that can be registered with the phase's iteratorFactoryMonitor.
	 * In this way this integrator is informed if the iteratorFactory object changes to another 
	 * instance in the phase.
	 */
	protected Observer iteratorFactoryObserver() {
	    return new Observer() {
	        //This is the action that is to happen if phase takes a new iteratorFactory
	        public void update(Observable o, Object arg) {
	            makeIterators((IteratorFactory)arg);
	        }
	    };
	}
  public synchronized void addIntervalListener(IntervalListener iil) {
    boolean added = false;
    //must check all possibilities because listener may implement multiple pbc interfaces
    if(iil instanceof IntervalListener.BeforePbc) {
        intervalListenersBeforePbc.addElement(iil);
        added = true;
    }
    if(iil instanceof IntervalListener.ImposePbc) {
        intervalListenersImposePbc.addElement(iil);
        added = true;
    }
    if(iil instanceof IntervalListener.AfterPbc) {
        intervalListenersAfterPbc.addElement(iil);
        added = true;
    }
    //if not implementing any of the pbc interfaces, default is afterPbc
    if(!added) intervalListenersAfterPbc.addElement(iil);
  }

  public synchronized void removeIntervalListener(IntervalListener iil) {
    boolean removed = false;
    if(iil instanceof IntervalListener.BeforePbc) {
        intervalListenersBeforePbc.removeElement(iil);
        removed = true;
    }
    if(iil instanceof IntervalListener.ImposePbc) {
        intervalListenersImposePbc.removeElement(iil);
        removed = true;
    }
    if(iil instanceof IntervalListener.AfterPbc) {
        intervalListenersAfterPbc.removeElement(iil);
        removed = true;
    }
    if(!removed) intervalListenersAfterPbc.removeElement(iil);
  }
  
    /**
     * Notifies registered listeners that an interval has passed.
     * Not synchronized, so unpredictable behavior if listeners are added
     * while notification is in process (this should be rare).
     */
  public void fireIntervalEvent(IntervalEvent iie) {
    iie.setBeforePbc(true);
    int n = intervalListenersBeforePbc.size();
    for(int i = 0; i < n; i++) {
        IntervalListener listener = (IntervalListener)intervalListenersBeforePbc.elementAt(i);
        listener.intervalAction(iie);
    }
    n = intervalListenersImposePbc.size();
    for(int i = 0; i < n; i++) {
        IntervalListener listener = (IntervalListener)intervalListenersImposePbc.elementAt(i);
        listener.intervalAction(iie);
    }
    iie.setBeforePbc(false);
    n = intervalListenersAfterPbc.size();
    for(int i = 0; i < n; i++) {
        IntervalListener listener = (IntervalListener)intervalListenersAfterPbc.elementAt(i);
        listener.intervalAction(iie);
    }
  }
  
  /**
   * Registers with the given integrator all listeners currently registered with
   * this integrator.  Removes all listeners from this integrator.
   */
  public synchronized void transferListenersTo(Integrator anotherIntegrator) {
    if(anotherIntegrator == this) return;
    int n = intervalListenersBeforePbc.size();
    for(int i = 0; i < n; i++) {
        IntervalListener listener = (IntervalListener)intervalListenersBeforePbc.elementAt(i);
        anotherIntegrator.addIntervalListener(listener);
    }
    n = intervalListenersImposePbc.size();
    for(int i = 0; i < n; i++) {
        IntervalListener listener = (IntervalListener)intervalListenersImposePbc.elementAt(i);
        anotherIntegrator.addIntervalListener(listener);
    }
    n = intervalListenersAfterPbc.size();
    for(int i = 0; i < n; i++) {
        IntervalListener listener = (IntervalListener)intervalListenersAfterPbc.elementAt(i);
        anotherIntegrator.addIntervalListener(listener);
    }
    intervalListenersBeforePbc.removeAllElements();
    intervalListenersImposePbc.removeAllElements();
    intervalListenersAfterPbc.removeAllElements();
  }
    
    /**
     * Accessor method for the number of integration steps to be
     * performed by this integrator after it is started.
     */
    public int getMaxSteps() {return maxSteps;}
    /**
     * Mutator method for the number of integration steps to be
     * performed by this integrator after it is started.  Can
     * be changed while integrator is running; if set to a value
     * less than number of steps already executed, integration will end.
     */
    public void setMaxSteps(int m) {maxSteps = m;}
    
    /**
     * Sets integrator to begin running on its own thread.  This is the normal
     * way to begin the integrator's activity.  Fires an event to listeners indicating
     * that integrator has started, calls the initialize method, and starts a new
     * thread that then enters the integrators run() method.  If integrator is already
     * running, method call return immediately and has no effect.
     */
    public void start() {
        if(running) return;
        running = true;
        fireIntervalEvent(new IntervalEvent(this, IntervalEvent.START));
        haltRequested = false;
//        isPaused = false;
        this.initialize();
        runner = new Thread(this);
        runner.start();
    }

    /**
     * Main loop for conduct of integration.  Repeatedly calls doStep() method,
     * while checking for halt/pause/reset requests, firing regular interval events,
     * and entering a brief sleep state if so indicated by doSleep flag.  Integration
     * loop continues until number of steps equals maxSteps field.  This method
     * is not the normal way to run the integrator, as it runs on the calling thread.
     * Instead, start should be used to set the integrator to run on its own new thread.
     */
    public void run() {
//        if(running) return; //do not allow two threads
        running = true;
        stepCount = 0;
        int iieCount = interval+1;
        while(stepCount < maxSteps) {
            while(pauseRequested) doWait();//keep this before resetRequest, since need for reset might naturally follow completion of pause
            if(resetRequested) {doReset(); resetRequested = false;}
            if(haltRequested) break;
            this.doStep();
            if(--iieCount == 0) {
                fireIntervalEvent(intervalEvent);
                iieCount = interval;
            }
            if(doSleep) {
                try { Thread.sleep(sleepPeriod); }
                catch (InterruptedException e) { }
            }
            stepCount++;
        }//end of while loop
        initialized = false;
        running = false;
        fireIntervalEvent(new IntervalEvent(this, IntervalEvent.DONE));
    }//end of run method
    
    /**
     * Requests that the integrator reset itself.  The actual action an
     * integrator takes to do this differs with the type of integrator.
     * The reset is not performed until the completion of the current 
     * integration step, or until the integrator is unpaused if currently
     * in a paused state.
     */
    public void reset() {resetRequested = true;}
    
    /**
     * Method to put integrator in a condition of being paused.
     */
    private synchronized void doWait() {
//        isPaused = true;
        notifyAll(); //release any threads waiting for pause to take effect
        try {
            wait(); //put in paused state
        } catch(InterruptedException e) {}
//        isPaused = false;
    }
    
    //suspend and resume functions
    /**
     * Requests that the integrator pause its execution.  The actual suspension
     * of execution occurs only after completion of the current integration step.
     * The calling thread is put in a wait state until the pause takes effect.
     */
    public synchronized void pause() {
        if(running && !pauseRequested/*!isPaused*/) {
            pauseRequested = true;
            try {
                wait();  //make thread requesting pause wait until pause is in effect
            } catch(InterruptedException e) {}
        }
    }
    /**
     * Removes the integrator from the paused state, resuming execution where it left off.
     */
    public synchronized void unPause() {pauseRequested = false; notifyAll();}
    /**
     * Queries whether the integrator is in a state of being paused.  This may
     * occur independent of whether the integrator is running or not.  If paused
     * but not running, then pause will take effect upon start.
     */
    public boolean isPaused() {return pauseRequested;}//isPaused;}
    
    /**
     * Indicates if the integrator has been started and has not yet completed.
     * If so, returns true, even if integrator is presently paused (but not halted).
     */
    public boolean isActive() {return running;}
    
    //stop function
    //consider having calling thread here join() or wait() for halt to take effect
    /**
     * Request that the integrator terminate its thread on the next integration step.
     * Does not cause calling thread to wait until this is completed, so it would
     * be prudent to have the calling thread join() to suspend it until the halt
     * is in effect.
     */
    public void halt() {
        if(running) haltRequested = true;
        if(pauseRequested) unPause();
    }
    
    /**
     * Method to make the calling thread wait until this integrator thread has finished.
     */
    public void join() {
        if(!running) return;
        try {
            runner.join();
        } catch(InterruptedException e) {}
    }
    
// Class generated by integrator as one of the properties of each atom

    interface Agent extends java.io.Serializable {
        /**
         * Integrator agent that holds a force vector.  Used to indicate that an atom could be
         * under the influence of a force.
         */
        interface Forcible extends Agent {
            public Space.Vector force();
        }
    }
    
    public static class IntervalEvent extends EventObject{
        
        //Typed constants used to indicate the type of event integrator is announcing
        public static final Type START = new Type("Start");       //simulation is starting
        public static final Type INTERVAL = new Type("Interval"); //routine interval event
        public static final Type DONE = new Type("Done");         //simulation is finished
        public static final Type INITIALIZE = new Type("Initialize"); //integrator is initializing
        
      //  private final Type type;//sometimes compiles, sometimes doesn't when declared final
        private Type type;
        private boolean beforePbc;
        
        public IntervalEvent(Integrator source, Type t) {
            super(source);
            type = t;
        }
        
        /**
         * Indicates if notification is before or after the ImposePbc listeners have been
         * notified.  Returns true if notifying BeforePbc or ImposePbc listeners;
         * return false if notifying AfterPbc listeners.
         */
        public final boolean isBeforePbc() {return beforePbc;}
        /**
         * Sets the before/after status relative to imposePbc listeners.
         * Should be used only by Integrator that is firing event.
         */
        final void setBeforePbc(boolean b) {beforePbc = b;}
        
        public Type type() {return type;}
        
        //class used to mark the different types of interval events
        private final static class Type extends Constants.TypedConstant {
           private Type(String label) {super(label);}
           public static final Constants.TypedConstant[] choices = new Constants.TypedConstant[] {START, INTERVAL, DONE, INITIALIZE};
           public final Constants.TypedConstant[] choices() {return choices;}
        }
    }
    
    public interface IntervalListener extends java.util.EventListener {
        public void intervalAction(IntervalEvent evt);
        
        /**
         * Marker interface that indicates an IntervalListener that should be notified before any
         * application of periodic boundary conditions.
         */
        public interface BeforePbc extends IntervalListener {}

        /**
         * Marker interface that indicates an IntervalListener that will invoke 
         * periodic boundary conditions when notified.
         */
        public interface ImposePbc extends IntervalListener {}

        /**
         * Marker interface that indicates an IntervalListener that should be notified after any
         * application of periodic boundary conditions.
         * This is the default for any IntervalListener that doesn't have a 
         * IntervalListenerPBC marker interface.
         */
        public interface AfterPbc extends IntervalListener {}
    }
    
}

