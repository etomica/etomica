package etomica;

/**
 * Parent class of all simulation elements, which are the classes of type
 * Species, Potential, Integrator, Controller, Phase, MeterAbstract, Device,
 * Display, or Action.  These are the classes that present the top-level 
 * programming interface to someone assembling a simulation applet or application.
 *
 * @author David Kofke
 */
public class SimulationElement implements java.io.Serializable {
    
    private final Simulation parentSimulation;
    public final int index;
    private boolean added = false;
    private final Class baseClass;
    String name;
    
    public SimulationElement(Simulation sim, Class baseClass) {
        parentSimulation = sim;
        this.baseClass = baseClass;
        index = sim.register(this);
        name = baseClass.getName().substring(8) + Integer.toString(index);
    }
    
    public final Simulation parentSimulation() {return parentSimulation;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
    public final Class baseClass() {return baseClass;}

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