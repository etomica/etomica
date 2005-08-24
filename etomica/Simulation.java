package etomica;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.action.activity.Controller;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.data.DataAccumulator;
import etomica.data.meter.Meter;
import etomica.integrator.Integrator;
import etomica.log.LoggerAbstract;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.util.Default;
import etomica.util.NameMaker;

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
 
public class Simulation extends EtomicaInfo implements java.io.Serializable  {
    
    /**
     * Class that implements the final tying up of the simulation elements before starting the simulation.
     * Default choice is the CoordinatorOneIntegrator.
     */
//    public Mediator elementCoordinator;
    protected HashMap elementLists = new HashMap(16);
    
    public final PotentialMaster potentialMaster;
    public final Space space;
    public final SpeciesRoot speciesRoot;
    
    
    /**
     * Constructs a default 2D, dynamic simulation.
     */
    public Simulation() {
        this(Space2D.getInstance());
    }
    
    /**
     * Creates a new simulation using the given space, with a default
     * setting of isDynamic = true.
     */
    public Simulation(Space space) {
        this(space, true, new PotentialMaster(space));
    }
    
    public Simulation(Space space, boolean isDynamic, PotentialMaster potentialMaster) {
        this(space, isDynamic, potentialMaster,Default.BIT_LENGTH);
    }
    
    public Simulation(Space space, boolean isDynamic, PotentialMaster potentialMaster, int[] bitLength) {
        this.space = space;
        this.dynamic = isDynamic;
        setName(NameMaker.makeName(this.getClass()));
//        elementCoordinator = new Mediator(this);
        this.potentialMaster = potentialMaster;
        setController(new Controller());
        speciesRoot = new SpeciesRoot((int[])bitLength.clone());
    }//end of constructor
                 
    public final Space space() {return space;}
    /**
     * @return the <code>nth</code> instantiated phase (indexing from zero)
     */
    public final Phase phase(int n) {return (Phase)getElement(Phase.class, n);}

    /**
     * @return the <code>nth</code> instantiated integrator (indexing from zero)
     */
    public final Integrator integrator(int n) {return (Integrator)getElement(Integrator.class, n);}
    /**
     * @return the <code>nth</code> instantiated meter (indexing from zero)
     */
    public final Meter meter(int n) {return (Meter)getElement(Meter.class, n);}
	/**
	 * @return the <code>nth</code> instantiated logger (indexing from zero)
	 */
	public final LoggerAbstract logger(int n) {return (LoggerAbstract)getElement(LoggerAbstract.class, n);}

	private Object getElement(Class clazz, int n) {
		LinkedList list = (LinkedList)elementLists.get(clazz);
		if(n < 0 || n >= list.size()) return null;
		return list.get(n);
	}
	
	public final LinkedList getLoggerList() {return loggerList;}
    public final LinkedList getDataAccumulatorList() {return dataAccumulatorList;}
    public final LinkedList getPhaseList() {
        phaseList.clear();
        AtomIteratorListSimple listIterator = new AtomIteratorListSimple(((AtomTreeNodeGroup)speciesRoot.node).childList);
        listIterator.reset();
        while(listIterator.hasNext()) {
            phaseList.add(listIterator.nextAtom().node.parentPhase());
        }
        return phaseList;
    }
    public final LinkedList getMeterList() {return meterList;}
    public final LinkedList getIntegratorList() {return integratorList;}
    public final LinkedList getSpeciesList() {return speciesList;}
    
    /**
     * Add the given DataAccumulator to a list kept by the simulation.
     * No other effect results from registering the DataAccumulator.  
     * The list of registered DataAccumulators may be retrieved via
     * the getDataManagerList method.  A DataAccumulator may be
     * removed from the list via the unregister method.
     */
    public void register(DataAccumulator dataManager) {
     	dataAccumulatorList.add(dataManager);
    }

    public void register(Phase phase) {
        phaseList.add(phase);
    }
    
    public void register(Meter meter) {
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
    public void unregister(DataAccumulator dataAccumulator) {
    	dataAccumulatorList.remove(dataAccumulator);
    }
     
    public void unregister(Phase phase) {
        phaseList.remove(phase);
    }
    
    public void unregister(Meter meter) {
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
    
    /**
     * @return Returns a flag indicating whether the simulation involves molecular dynamics.
     */
    public boolean isDynamic() {
        return dynamic;
    }
    
    public static final java.util.Random random = new java.util.Random();
//    public static final java.util.Random random = new java.util.Random(1);
        
     private final boolean dynamic;
     private Controller controller;     
     private final LinkedList dataAccumulatorList = new LinkedList();
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


