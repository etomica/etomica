package etomica;

import etomica.units.UnitSystem;
import etomica.log.LoggerAbstract;

//Java2 imports
//import java.util.HashMap;
//import java.util.LinkedList;
//import java.util.Iterator;

import etomica.utility.java2.HashMap;
import etomica.utility.java2.LinkedList;
import etomica.utility.java2.Iterator;

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
 
public class Simulation extends SimulationElement {
    
    /**
     * Flag indicating whether simulation is being run within Etomica editor application.
     * This is set to true by Etomica if it is running; otherwise it is false.
     */
    public static boolean inEtomica = false;
    /**
     * Class that implements the final tying up of the simulation elements before starting the simulation.
     * Default choice is the CoordinatorOneIntegrator.
     */
    public Mediator elementCoordinator;
    protected HashMap elementLists = new HashMap(16);
    
    public final PotentialMaster potentialMaster;
    
    private static final LinkedList instances = new LinkedList();
    
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
        super(space, instanceCount++, Simulation.class);
        instance = this;
        instances.add(this);
        setName("Simulation" + Integer.toString(instanceCount));
        elementLists.put(Species.class, new LinkedList());
        elementLists.put(Integrator.class, new LinkedList());
        elementLists.put(Phase.class, new LinkedList());
        elementLists.put(Controller.class, new LinkedList());
		elementLists.put(MeterAbstract.class, new LinkedList());
		elementLists.put(LoggerAbstract.class, new LinkedList());
        elementCoordinator = new Mediator(this);
        this.potentialMaster = potentialMaster;
        instantiationEventManager.fireEvent(new SimulationEvent(this));
    }//end of constructor
                 
    public final Space space() {return space;}
    /**
     * @return the <code>nth</code> instantiated phase (indexing from zero)
     */
    public final Phase phase(int n) {return (Phase)phaseList().get(n);}
    /**
     * @return the <code>nth</code> instantiated species (indexing from zero)
     */
    public final Species species(int n) {return (Species)speciesList().get(n);}
    /**
     * @return the <code>nth</code> instantiated controller (indexing from zero)
     */
    public final Controller controller(int n) {return (Controller)controllerList().get(n);}
    /**
     * @return the <code>nth</code> instantiated integrator (indexing from zero)
     */
    public final Integrator integrator(int n) {return (Integrator)integratorList().get(n);}
    /**
     * @return the <code>nth</code> instantiated meter (indexing from zero)
     */
    public final MeterAbstract meter(int n) {return (MeterAbstract)meterList().get(n);}
	/**
	 * @return the <code>nth</code> instantiated logger (indexing from zero)
	 */
	public final LoggerAbstract logger(int n) {return (LoggerAbstract)loggerList().get(n);}

    public final LinkedList phaseList() {return (LinkedList)elementLists.get(Phase.class);}
	public final LinkedList meterList() {return (LinkedList)elementLists.get(MeterAbstract.class);}
	public final LinkedList loggerList() {return (LinkedList)elementLists.get(LoggerAbstract.class);}
    public final LinkedList speciesList() {return (LinkedList)elementLists.get(Species.class);}
    public final LinkedList integratorList() {return (LinkedList)elementLists.get(Integrator.class);}
    public final LinkedList controllerList() {return (LinkedList)elementLists.get(Controller.class);}
    public final LinkedList elementList(Class clazz) {return (LinkedList)elementLists.get(clazz);}
    
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
    
    int register(SimulationElement element) {
    	super.register(element);
        LinkedList list = (LinkedList)elementLists.get(element.baseClass());
        if(list.contains(element)) return -1;
        list.add(element);
        return list.size() - 1;
    }//end of register method
                
    public void unregister(SimulationElement element) {
    	super.unregister(element);
        LinkedList list = (LinkedList)elementLists.get(element.baseClass());
        if(!list.contains(element)) return;
        list.remove(element);
    }
     
    public void resetIntegrators() {
        for(Iterator is=integratorList().iterator(); is.hasNext(); ) {
            Integrator integrator = (Integrator)is.next();
            integrator.reset();
        }
    }
    
    private static int instanceCount = 0;
    
//    public static final java.util.Random random = new java.util.Random();
    public static final java.util.Random random = new java.util.Random(1);
        
     /**
      * Returns the mediator that coordinates the elements of the simulation.
      * The is the same as the elementCoordinator field, but provides another
      * way to access it.  This method may someday supercede direct access to
      * the elementCoordinator field, so it is the preferred way to access it.
      */
     public Mediator mediator() {return elementCoordinator;}
     
     public Simulation simulation() {return this;}
     
     public static final SimulationEventManager instantiationEventManager = new SimulationEventManager();
 
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
	    Controller controller = new Controller();
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


