//This class includes a main method to demonstrate its use
package etomica;

import etomica.units.UnitSystem;

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
public class Simulation extends SimulationElement {
    
    public static final String VERSION = "Simulation:01.11.20";
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
    
    public final Hamiltonian hamiltonian;
    
    public PotentialCalculationEnergySum energySum;
    
   /**
    * Object describing the nature of the physical space in which the simulation is performed
    */
 //   public final Space space;
    

    //default unit system for I/O (internal calculations are all done in simulation units)
    private static UnitSystem unitSystem = new UnitSystem.Sim();
    
    /**
     * A static instance of a Simulation, for which the current value at any time is
     * used as a default simulation in many places.  Any new instance of a Simulation
     * is assigned to this field upon construction.
     */
    public static Simulation instance = new Simulation(new Space2D());
       
    public Simulation() {
        this(new Space2D());
    }
    
    /**
     * Creates a new simulation using the given space, and sets the static
     * instance field equal to the new instance.
     */
    public Simulation(Space s) {
        super(s, instanceCount++, Simulation.class);
        instance = this;
        setName("Simulation" + Integer.toString(instanceCount++));
        elementLists.put(Potential.class, new LinkedList());
        elementLists.put(Species.class, new LinkedList());
        elementLists.put(Integrator.class, new LinkedList());
        elementLists.put(Phase.class, new LinkedList());
        elementLists.put(Controller.class, new LinkedList());
        elementLists.put(MeterAbstract.class, new LinkedList());
        elementCoordinator = new Mediator(this);
        hamiltonian = new Hamiltonian(this);
    }//end of constructor
    
    /**
     * Accessor method for the default I/O unit system.
     */
    public static final UnitSystem unitSystem() {return unitSystem;}
    /**
     * Accessor method for the default I/O unit system.
     */
    public static final void setUnitSystem(UnitSystem us) {unitSystem = us;}
              
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
     * @return the <code>nth</code> instantiated potential (indexing from zero)
     */
    public final Potential potential(int n) {return (Potential)potentialList().get(n);}
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

    public final LinkedList phaseList() {return (LinkedList)elementLists.get(Phase.class);}
    public final LinkedList meterList() {return (LinkedList)elementLists.get(MeterAbstract.class);}
    public final LinkedList speciesList() {return (LinkedList)elementLists.get(Species.class);}
    public final LinkedList integratorList() {return (LinkedList)elementLists.get(Integrator.class);}
    public final LinkedList controllerList() {return (LinkedList)elementLists.get(Controller.class);}
    public final LinkedList potentialList() {return (LinkedList)elementLists.get(Potential.class);}
    public final LinkedList elementList(Class clazz) {return (LinkedList)elementLists.get(clazz);}
  
    int register(SimulationElement element) {
    	super.register(element);
//        if(hamiltonian == null || element == hamiltonian.potential) return;
        LinkedList list = (LinkedList)elementLists.get(element.baseClass());
        if(list.contains(element)) return -1;
//        if(element instanceof Potential && !(element instanceof Potential.Null))
//                hamiltonian.potential.addPotential((Potential)element);
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
    
    /**
     * Returns the energy-sum object used by all Integrator, Meters, MCMoves, etc.
     * Returns the object last passed to setEnergySum, or if no such call was
     * made (or if the call set energySum to null) a new instance of 
     * PotentialCalculationEnergySum is returned with each call.
     */
    public PotentialCalculationEnergySum energySum(Object obj) {
        if(energySum != null) return energySum;
        else return new PotentialCalculationEnergySum();
    }
    /**
     * Sets the energy-sum object that will be given to all Integrators, Meters, etc.
     * that require one for their operation.  Other such instances may perform additional
     * calculations with the energy measurement. May be set to null.
     */
    public void setEnergySum(PotentialCalculationEnergySum sum) {
        energySum = sum;
    }
    
    private static int instanceCount = 0;
    
    public static final java.util.Random random = new java.util.Random();
//    public static final java.util.Random random = new java.util.Random(1);
        
     /**
      * Returns the mediator that coordinates the elements of the simulation.
      * The is the same as the elementCoordinator field, but provides another
      * way to access it.  This method may someday supercede direct access to
      * the elementCoordinator field, so it is the preferred way to access it.
      */
     public Mediator mediator() {return elementCoordinator;}
     
     public void setIteratorFactory(IteratorFactory factory) {
        iteratorFactory = factory;
     }
     public IteratorFactory getIteratorFactory() {return iteratorFactory;}
     
     public IteratorFactory iteratorFactory = IteratorFactorySimple.INSTANCE;
     
     public Simulation simulation() {return this;}
 
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


