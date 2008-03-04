package etomica.simulation;

import etomica.api.IBox;
import etomica.api.IController;
import etomica.api.IRandom;
import etomica.api.ISimulation;

import etomica.action.activity.Controller;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.util.Arrays;
import etomica.util.RandomNumberGenerator;

/**
 * The main class that organizes the elements of a molecular simulation.
 * Holds a single Space instance that is referenced in
 * many places to obtain spatial elements such as vectors and boundaries.
 */

public class Simulation implements java.io.Serializable, ISimulation  {


    /**
     * Creates a new simulation using the given space, with a default
     * setting of isDynamic = true.
     */
    public Simulation(Space space) {
        this(space, true);
    }
    
    public Simulation(Space space, boolean isDynamic) {
        this.space = space;
        this.dynamic = isDynamic;
        boxList = new IBox[0];
        setController(new Controller());
        random = new RandomNumberGenerator();
        eventManager = new SimulationEventManager();
        speciesManager = new SpeciesManager(this);
    }

    public final void addBox(IBox newBox) {
        for (int i=0; i<boxList.length; i++) {
            if (boxList[i] == newBox) {
                throw new IllegalArgumentException("Box "+newBox+" is already a part of this Simulation");
            }
        }
        boxList = (IBox[])Arrays.addObject(boxList, newBox);
        newBox.resetIndex(this);
        speciesManager.boxAddedNotify(newBox);
        eventManager.fireEvent(new SimulationBoxAddedEvent(newBox));
    }
    
    public final void removeBox(IBox oldBox) {
        boolean found = false;
        for (int i=0; i<boxList.length; i++) {
            if (boxList[i] == oldBox) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw new IllegalArgumentException("Box "+oldBox+" is not part of this Simulation");
        }

        boxList = (IBox[])Arrays.removeObject(boxList, oldBox);

        for (int i = oldBox.getIndex(); i<boxList.length; i++) {
            boxList[i].resetIndex(this);
        }
        
        // notify oldBox that we no longer have it.
        oldBox.resetIndex(null);
        
        eventManager.fireEvent(new SimulationBoxRemovedEvent(oldBox));
    }
    
    /**
     * Returns an array of Boxs contained in the Simulation
     */
    public final IBox[] getBoxs() {
        return boxList;
    }

    /**
     * Returns the Controller used to run the simulation's Actions and 
     * Activities.
     */
	public IController getController() {
		return controller;
	}
	
	//TODO transfer control from old to new controller (copy over integrators, etc)
    //AJS really?
	public void setController(IController controller) {
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
    private IBox[] boxList;
    private final SpeciesManager speciesManager;
    protected final IRandom random;
    protected final boolean dynamic;
    private IController controller;
}
