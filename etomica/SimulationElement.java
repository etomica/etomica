package etomica;
import etomica.utility.java2.LinkedList;


/**
 * Parent class of all simulation elements, which are the classes of type
 * Simulation, Species, Potential, Integrator, Controller, Phase, MeterAbstract,
 * Device, Display, or Action.  These are the classes that present the top-level
 * programming interface to someone assembling a simulation applet or
 * application.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 07/03/02 (DAK/SKK) Added reset method 
  * 01/27/03 (DAK) Revised to act as
  * container for other simulation elements, and to be superclass of Simulation
  * 
  */
  
public class SimulationElement implements java.io.Serializable {
    
    public final SimulationElement parentElement;
    public final Space space;
    public final int index;
    private boolean added = false;
    private final Class baseClass;
    protected String name;
    
    /**
     * Constructor for use only by a Simulation subclass.
     * @param space
     * @param index
     * @param baseClass
     */
    SimulationElement(Space space, int index, Class baseClass) {
    	if(!(this instanceof Simulation)) throw new RuntimeException("Inappropriate use of SimulationElement constructor");
    	this.space = space;
    	this.index = index;
    	this.baseClass = baseClass;
    	parentElement = null;
    }
    protected SimulationElement(SimulationElement parentElement, Class baseClass) {
        this.parentElement = parentElement;
        space = parentElement.space;
        this.baseClass = baseClass;
        index = parentElement.register(this);
        name = baseClass.getName().substring(8) + Integer.toString(index);
    }
        
    public SimulationElement parentElement() {return parentElement;}
    public Simulation simulation() {return parentElement.simulation();}
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
     * Method to set the name of this simulation element. The element's name
     * provides a convenient way to label output data that is associated with
     * it.  This method might be used, for example, to place a heading on a
     * column of data. Default name is the base class followed by the integer
     * index of this element.
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
    
    int register(SimulationElement element) {
    	if(element.parentElement != this) throw new IllegalArgumentException("Cannot invoke ");
		//add to list of all elements.  Put Species and Phase classes first in list
		//so that all species agents are deployed in phases by mediator before it processes
		//other elements.
		if(element.baseClass() == Phase.class || 
			element.baseClass() == Species.class) allElements.addFirst(element);
		else allElements.addLast(element);
		return allElements.size() - 1;
    }
    
	public void unregister(SimulationElement element) {
		allElements.remove(element);
	}
 
	/**
	 * Returns a clone of the list of all elements added to this element.
	 * @return LinkedList  Copy of the current list of all elements.
	 */
	public LinkedList allElements() {
		LinkedList list;
 //       synchronized(this) {
			list = (LinkedList)allElements.clone();
 //       }
		return list;
	}
    
	/**
	 * List of all simulation elements contained by this one.
	 */
	 private final LinkedList allElements = new LinkedList();
          
}