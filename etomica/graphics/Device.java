package etomica.graphics;
import java.awt.Component;

import etomica.Action;
import etomica.Controller;
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
    protected Controller controller;
    
    public Device() {
        this(null);
    }
    
    public Device(Controller controller) {
        this.controller = controller;
        setName(NameMaker.makeName(this.getClass()));
    }
    
    /**
     * Method called by subclasses to execute the action of the device.
     * Action is specified by the subclass as the targetAction field; if
     * this field is null, not action is performed.
     * The action is normally performed by calling the doActionNow method
     * of the controller, which causes the controller's thread to perform
     * the action as soon as possible after suspending all simulation activities.
     * If a controller has not been defined (is null), then the action is performed
     * on the thread calling this method.
     */
    protected void doAction(Action action) {
        if (action == null) return;
        if(controller != null) {
            controller.doActionNow(action);
        } else {
            action.actionPerformed();
        }
    }
      
    
    /**
     * @return Returns the controller.
     */
    public Controller getController() {
        return controller;
    }
    /**
     * @param controller The controller to set.
     */
    public void setController(Controller controller) {
        this.controller = controller;
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