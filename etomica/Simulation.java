package etomica;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.log.LoggerAbstract;
import etomica.space2d.Space2D;
import etomica.utility.NameMaker;

/**
 * The main class that organizes the elements of a molecular simulation.
 * Holds a single space object that is referenced in
 * many places to obtain spatial elements such as vectors and boundaries.  Also
 * holds an object that specifies the unit system used to default all I/O.  A single
 * instance of Simulation is held as a static field, and which forms the default
 * Simulation class needed in the constructor of all simulation elements.
 *
 * @author David Kofke
 */

/* History
 * 04/18/03 (DAK) added LoggerList
 */
 
public class Simulation {
    
    /**
     * Flag indicating whether simulation is being run within Etomica editor application.
     * This is set to true by Etomica if it is running; otherwise it is false.
     */
    public static boolean inEtomica = false;
    /**
     * Class that implements the final tying up of the simulation elements before starting the simulation.
     * Default choice is the CoordinatorOneIntegrator.
     */
//    public Mediator elementCoordinator;
    protected HashMap elementLists = new HashMap(16);
    
    public final PotentialMaster potentialMaster;
    public final Space space;
    
    /**
     * A static instance of a Simulation, for which the current value at any time is
     * used as a default simulation in many places.  Any new instance of a Simulation
     * is assigned to this field upon construction.
     */
    public static Simulation instance;// = new Simulation(new Space2D());
       
    public Simulation() {
        this(new Space2D());
    }
    
    /**
     * Creates a new simulation using the given space, and sets the static
     * instance field equal to the new instance.
     */
    public Simulation(Space space) {
    	this(space, new PotentialMaster(space));
    }
    
    public Simulation(Space space, PotentialMaster potentialMaster) {
        this.space = space;
        instance = this;
        instanceCount++;
        instances.add(this);
        setName(NameMaker.makeName(this.getClass()));
//        elementCoordinator = new Mediator(this);
        this.potentialMaster = potentialMaster;
        potentialMaster.setSimulation(this);
        setController(new Controller());
        instantiationEventManager.fireEvent(new SimulationEvent(this));
    }//end of constructor
                 
    public final Space space() {return space;}
    /**
     * @return the <code>nth</code> instantiated phase (indexing from zero)
     */
    public final Phase phase(int n) {return (Phase)getElement(Phase.class, n);}
    /**
     * @return the <code>nth</code> instantiated species (indexing from zero)
     */
    public final Species species(int n) {return (Species)getElement(Species.class, n);}
    /**
     * @return the <code>nth</code> instantiated controller (indexing from zero)
     * @deprecated use getController instead
     */
    public final Controller controller(int n) {
    	if(n != 0) throw new IllegalArgumentException("only one controller; must call with 0 index");
    	return getController();
    }
    /**
     * @return the <code>nth</code> instantiated integrator (indexing from zero)
     */
    public final Integrator integrator(int n) {return (Integrator)getElement(Integrator.class, n);}
    /**
     * @return the <code>nth</code> instantiated meter (indexing from zero)
     */
    public final MeterAbstract meter(int n) {return (MeterAbstract)getElement(MeterAbstract.class, n);}
	/**
	 * @return the <code>nth</code> instantiated logger (indexing from zero)
	 */
	public final LoggerAbstract logger(int n) {return (LoggerAbstract)getElement(LoggerAbstract.class, n);}

	private Object getElement(Class clazz, int n) {
		LinkedList list = (LinkedList)elementLists.get(clazz);
		if(n < 0 || n >= list.size()) return null;
		else return list.get(n);
	}
	
	public final LinkedList getLoggerList() {return loggerList;}
    public final LinkedList getDataManagerList() {return dataManagerList;}
    public final LinkedList getPhaseList() {return phaseList;}
    public final LinkedList getMeterList() {return meterList;}
    public final LinkedList getIntegratorList() {return integratorList;}
    public final LinkedList getSpeciesList() {return speciesList;}
    
    public static Simulation getDefault() {
    	if(instance == null) instance = new Simulation(new Space2D());
    	return instance;
    }
    
    public static void setDefault(Simulation sim) {
    	instance = sim;
    }
    
    public static LinkedList getInstances() {
    	return instances;
    }
  
    public void dispose() {
    	instances.remove(this);
    }
    
    /**
     * Add the given dataManger to a list kept by the simulation.
     * No other effect results from registering the dataManager.  
     * The list of registered accumulationManagers may be retrieved via
     * the getDataManagerList method.  An dataManager may be
     * removed from the list via the unregister method.
     * @param dataManager
     */
    public void register(DataManager dataManager) {
     	dataManagerList.add(dataManager);
    }

    public void register(Phase phase) {
        phaseList.add(phase);
    }
    
    public void register(MeterAbstract meter) {
        meterList.add(meter);
    }
    
    public void register(Integrator integrator) {
        integratorList.add(integrator);
    }
    
    public void register(Species species) {
        speciesList.add(species);
    }
    
    public void register(LoggerAbstract logger) {
        loggerList.add(logger);
    }
    
    /**
     * Removes the given dataManager from the list of dataManagers
     * kept by the simulation.  No other action results upon removing it from
     * this list.  If the given dataManager is not in the list already,
     * the method returns without taking any action.
     */
    public void unregister(DataManager dataManager) {
    	dataManagerList.remove(dataManager);
    }
     
    public void unregister(Phase phase) {
        phaseList.remove(phase);
    }
    
    public void unregister(MeterAbstract meter) {
        meterList.remove(meter);
    }
    
    public void unregister(Integrator integrator) {
        integratorList.remove(integrator);
    }
    
    public void unregister(Species species) {
        speciesList.remove(species);
    }
    
    public void unregister(LoggerAbstract logger) {
        loggerList.remove(logger);
    }
    
    public void resetIntegrators() {
        for(Iterator is=getIntegratorList().iterator(); is.hasNext(); ) {
            Integrator integrator = (Integrator)is.next();
            integrator.reset();
        }
    }
    
	public Controller getController() {
		return controller;
	}
	
	//TODO transfer control from old to new controller (copy over integrators, etc)
	public void setController(Controller controller) {
		this.controller = controller;
	}
    
    /**
     * Accessor method of the name of this simulation.
     * 
     * @return The given name of this phase
     */
    public final String getName() {return name;}
    /**
     * Method to set the name of this simulation. The simulation's name
     * provides a convenient way to label output data that is associated with
     * it.  This method might be used, for example, to place a heading on a
     * column of data. Default name is the base class followed by the integer
     * index of this simulation.
     * 
     * @param name The name string to be associated with this element
     */
    public void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the phase
     */
    public String toString() {return getName();}
    
    private static int instanceCount = 0;
    
    public static final java.util.Random random = new java.util.Random();
//    public static final java.util.Random random = new java.util.Random(1);
        
     /**
      * Returns the mediator that coordinates the elements of the simulation.
      * The is the same as the elementCoordinator field, but provides another
      * way to access it.  This method may someday supercede direct access to
      * the elementCoordinator field, so it is the preferred way to access it.
      */
//     public Mediator mediator() {return elementCoordinator;}
     
     public Simulation simulation() {return this;}
     
     public static final SimulationEventManager instantiationEventManager = new SimulationEventManager();
     private static final LinkedList instances = new LinkedList();
     private Controller controller;     
     private final LinkedList dataManagerList = new LinkedList();
     private final LinkedList phaseList = new LinkedList();
     private final LinkedList loggerList = new LinkedList();
     private final LinkedList meterList = new LinkedList();
     private final LinkedList integratorList = new LinkedList();
     private final LinkedList speciesList = new LinkedList();
     private String name;
     /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
        Default.ATOM_SIZE = 1.0;                   
	    IntegratorHard integratorHard = new IntegratorHard();
	    SpeciesSpheres speciesSpheres = new SpeciesSpheres();
	    speciesSpheres.setNMolecules(300);
	    Phase phase = new Phase();
	    Potential2 potential = new P2HardSphere();
	    DisplayPhase displayPhase = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard.new Timer(integratorHard.chronoMeter());
	    timer.setUpdateInterval(10);
        integratorHard.setTimeStep(0.01);
        displayPhase.setColorScheme(new ColorSchemeNull());
        for(Atom atom=phase.firstAtom(); atom!=null; atom=atom.nextAtom()) {
            atom.setColor(Constants.randomColor());
        }
        
        //this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		                                    
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);
        Simulation.makeAndDisplayFrame(Simulation.instance);
        
     //   controller.start();
    }//end of main
 */   
}


