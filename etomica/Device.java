package etomica;
import etomica.units.Unit;
import etomica.units.Dimensioned;
import etomica.units.Dimension;
import java.awt.Component;

/**
 * Base class for all Devices.  These are objects that permit manipulation of the
 * fields of other objects.  They provide a (usually) graphical interface to permit
 * values of fields to be input and changed.
 */
public abstract class Device implements Simulation.GraphicalElement, Dimensioned, java.io.Serializable {
    
    public static final String VERSION = "Device:01.01.17";
    
    protected Unit unit;
    private String name;
    private final Simulation parentSimulation;
    private boolean added = false;
    
    public Device(Simulation sim) {
        parentSimulation = sim;
        parentSimulation.register(this);
    }
    
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Device.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
    public abstract Component graphic(Object obj);
    
    public void setUnit(Unit u) {unit = u;}
    public Unit getUnit() {return unit;}
    public Dimension dimension() {return unit.dimension();} //may want to override this in most cases

    /**
     * Accessor method of the name of this object
     * 
     * @return The given name
     */
    public final String getName() {return name;}

    /**
     * Method to set the name of this object
     * 
     * @param name The name string to be associated with this object
     */
    public final void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the object
     */
    public String toString() {return getName();}  //override Object method
          
    
}