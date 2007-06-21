package etomica.simulation;

import etomica.action.activity.Controller;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.util.Arrays;
import etomica.util.IRandom;
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

    private static final long serialVersionUID = 4L;
    protected final Space space;
    protected final SimulationEventManager eventManager;
    private Phase[] phaseList;
    private final SpeciesManager speciesManager;
    protected final IRandom random;
    protected final boolean dynamic;
    private Controller controller;
}
