package etomica;

/**
 * Parent class of all simulation elements, which are the classes of type
 * Species, Potential, Integrator, Controller, Phase, MeterAbstract, Device,
 * Display, or Action.  These are the classes that present the top-level 
 * programming interface to someone assembling a simulation applet or application.
 *
 * @author David Kofke
 */
 
 /* History of changes
  *  7/3/02 (DAK/SKK) Added reset method 
  */
  
public class SimulationElement implements java.io.Serializable {
    
    private final SimulationElement parentElement;
    public final Space space;
    public final int index;
    private boolean added = false;
    private final Class baseClass;
    protected String name;
    
    public SimulationElement(Simulation sim, Class baseClass) {
        parentSimulation = sim;
        space = sim.space;
        this.baseClass = baseClass;
        index = sim.register(this);
        name = baseClass.getName().substring(8) + Integer.toString(index);
    }
    
    /**
     * Constructor for situtions when element is not to be used
     * directly as part of a simulation.
     */
    public SimulationElement(Space space, Class baseClass, int index) {
        parentSimulation = null;
        this.space = space;
        this.baseClass = baseClass;
        this.index = index;
    }
    
    public final Simulation parentSimulation() {return parentSimulation;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
    public final Class baseClass() {return baseClass;}
    
    /**
     * Resets the element to some initial condition.  Invoking this on all elements in
     * a simulation should put the system in an "initialized" state, ready to begin
     * a new simulation run.  It is not guaranteed that the initial state will be identical
     * to that when the simulation is loaded.  Changes in the instance of boundary, system size,
     * the instance of configuration in phase, and so on, are not in general undone to restore
     * the system to its original condition.
     */
    public void reset() {}

    /**
     * Accessor method of the name of this species
     * 
     * @return The given name of this species
     */
    public final String getName() {return name;}
    /**
     * Method to set the name of this species
     * The element's name provides a convenient way to label output data that is 
     * associated with it.  This method might be used, for example, to place
     * a heading on a column of data.
     * Default name is the base class followed by the integer index of this element.
     * 
     * @param name The name string to be associated with this element
     */
    public void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the species
     */
    public String toString() {return getName();}
          
}