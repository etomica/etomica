package etomica.graphics;
import java.awt.Component;

import etomica.units.Dimension;
import etomica.units.Dimensioned;
import etomica.units.Unit;
import etomica.utility.NameMaker;

/**
 * Base class for all Devices.  These are objects that permit manipulation of the
 * fields of other objects.  They provide a (usually) graphical interface to permit
 * values of fields to be input and changed.
 */
public abstract class Device implements GraphicalElement, Dimensioned, java.io.Serializable {
    
    protected Unit unit;
    
    public Device() {
        setName(NameMaker.makeName(this.getClass()));
    }
        
    public abstract Component graphic(Object obj);
    
    /**
     * Same as graphic method with a null argument.
     */
    public final Component graphic() {return graphic(null);}
    
    public void setUnit(Unit u) {unit = u;}
    public Unit getUnit() {return unit;}
    public Dimension dimension() {return unit.dimension();} //may want to override this in most cases
    /**
     * Accessor method of the name of this device
     * 
     * @return The given name of this device
     */
    public final String getName() {return name;}
    /**
     * Method to set the name of this device
     */
    public void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the device
     */
    public String toString() {return getName();}

    private String name;
}