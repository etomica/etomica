package etomica.simulation;

import java.util.LinkedList;

import etomica.action.activity.Controller;
import etomica.data.DataSource;
import etomica.integrator.IIntegrator;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.util.Arrays;
import etomica.util.IRandom;
import etomica.util.NameMaker;
import etomica.util.RandomNumberGenerator;

/**
 * The main class that organizes the elements of a molecular simulation.
 * Holds a single Space instance that is referenced in
 * many places to obtain spatial elements such as vectors and boundaries.
 */

public class Simulation implements java.io.Serializable, ISimulation  {

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
        this(space, true);
    }
    
    public Simulation(Space space, boolean isDynamic) {
        this(space, isDynamic, new int[] {4, 19, 9});
    }
    
    public Simulation(Space space, boolean isDynamic, int[] bitLength) {
        this.space = space;
        this.dynamic = isDynamic;
        phaseList = new Phase[0];
        setName(NameMaker.makeName(this.getClass()));
        setController(new Controller());
        random = new RandomNumberGenerator();
        eventManager = new SimulationEventManager();
        speciesManager = new SpeciesManager(this, bitLength);
    }

    public final void addPhase(Phase newPhase) {
        for (int i=0; i<phaseList.length; i++) {
            if (phaseList[i] == newPhase) {
                throw new IllegalArgumentException("Phase "+newPhase+" is already a part of this Simulation");
            }
        }
        phaseList = (Phase[])Arrays.addObject(phaseList, newPhase);
        newPhase.resetIndex(this);
        speciesManager.phaseAddedNotify(newPhase);
        eventManager.fireEvent(new SimulationPhaseAddedEvent(newPhase));
    }
    
    public final void removePhase(Phase oldPhase) {
        boolean found = false;
        for (int i=0; i<phaseList.length; i++) {
            if (phaseList[i] == oldPhase) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw new IllegalArgumentException("Phase "+oldPhase+" is not part of this Simulation");
        }

        phaseList = (Phase[])Arrays.removeObject(phaseList, oldPhase);

        for (int i = oldPhase.getIndex(); i<phaseList.length; i++) {
            phaseList[i].resetIndex(this);
        }
        
        // notify oldPhase that we no longer have it.
        oldPhase.resetIndex(null);
        
        eventManager.fireEvent(new SimulationPhaseRemovedEvent(oldPhase));
    }
    
    /**
     * Returns an array of Phases contained in the Simulation
     */
    public final Phase[] getPhases() {
        return phaseList;
    }

    /**
     * Returns a list of Integrators that have been registered with the Simulation
     */
    public final LinkedList getIntegratorList() {return integratorList;}
    
    /**
     * Returns an array of DataStreamsHeaders for data streams that have been
     * registered with the simulation.
     */
    public DataStreamHeader[] getDataStreams() {
        return (DataStreamHeader[])dataStreams.clone();
    }
    
    /**
     * Clears the list of data streams that have been registered with the 
     * Simulation.
     */
    public void clearDataStreams() {
        dataStreams = new DataStreamHeader[0];
    }
    
    /**
     * Clears the list of Integrators that have been registered with the 
     * Simulation.
     */
    public void clearIntegrators() {
        integratorList.clear();
    }
    
    /**
     * Add the given Integrator to a list kept by the simulation.  No other 
     * effect results from registering the Integrator. The list of registered 
     * Integrators may be retrieved via the getIntegrators method.  An 
     * individual Integrator may be removed from the list via the unregister 
     * method.
     */
    public void register(IIntegrator integrator) {
        integratorList.add(integrator);
    }

    /**
     * Add the given DataSource and client (the object which calls the 
     * DataSource's getData method) to a list of data streams kept by the 
     * simulation.  No other effect results from registering the data stream. 
     * The list of registered data streams may be retrieved via the 
     * getDataStreams method.  A single data stream may be removed from the 
     * list via the unregister method.
     */
    public void register(DataSource dataSource, Object client) {
        for (int i=0; i<dataStreams.length; i++) {
            if (dataStreams[i].getDataSource() == dataSource && dataStreams[i].getClient() == client) {
                return;
            }
        }
        dataStreams = (DataStreamHeader[])Arrays.addObject(dataStreams,new DataStreamHeader(dataSource,client));
    }
    
    /**
     * Removes the given Integrator from the list of Integrators kept by the 
     * simulation.  No other action results upon removing it from this list.  
     * If the given Integrator is not in the list already, the method returns 
     * without taking any action.
     */
    public void unregister(IIntegrator integrator) {
        integratorList.remove(integrator);
    }
    
    /**
     * Removes the given data stream that represents the given DataSource and 
     * client from the list of data streams kept by the simulation.  No other 
     * action results upon removing it from this list.  If the given data 
     * stream has not been registered, the method returns without taking any 
     * action.
     */
    public void unregister(DataSource dataSource, Object client) {
        for (int i=0; i<dataStreams.length; i++) {
            if (dataStreams[i].getDataSource() == dataSource && dataStreams[i].getClient() == client) {
                dataStreams = (DataStreamHeader[])Arrays.removeObject(dataStreams,dataStreams[i]);
                return;
            }
        }
    }
    
    /**
     * Returns the Controller used to run the simulation's Actions and 
     * Activities.
     */
	public Controller getController() {
		return controller;
	}
	
	//TODO transfer control from old to new controller (copy over integrators, etc)
    //AJS really?
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
    
    /**
     * @return the space
     */
    public final Space getSpace() {
        return space;
    }

    public IRandom getRandom() {
        return random;
    }
    
    public SimulationEventManager getEventManager() {
        return eventManager;
    }
    
    public SpeciesManager getSpeciesManager() {
        return speciesManager;
    }

    private static final long serialVersionUID = 3L;
    protected final Space space;
    protected final SimulationEventManager eventManager;
    private Phase[] phaseList;
    private final SpeciesManager speciesManager;
    protected final IRandom random;
    private final boolean dynamic;
    private Controller controller;     
    private final LinkedList integratorList = new LinkedList();
    private DataStreamHeader[] dataStreams = new DataStreamHeader[0];
    private String name;
}
