package etomica.graphics;
import java.awt.Component;

import etomica.Controller;
import etomica.action.Action;
import etomica.units.Dimension;
import etomica.units.Dimensioned;
import etomica.units.Unit;
import etomica.util.NameMaker;

/**
 * Base class for all Devices.  These are objects that permit manipulation of the
 * fields of other objects.  They provide a (usually) graphical interface to permit
 * values of fields to be input and changed.
 */
public abstract class Device implements GraphicalElement, Dimensioned, java.io.Serializable {
    
    protected Unit unit;
    protected Controller controller;
    private final ActionSet actionSet = new ActionSet();
    
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
        actionSet.action = action;
        if(controller != null) {
            controller.doActionNow(actionSet);
        } else {
            actionSet.actionPerformed();
        }
    }
    
    /**
     * Sets an action to be performed before the primary action executed
     * by the device.  This is useful for modifying the behavior of concrete
     * subclasses of the device without accessing their primary action.
     */
    public void setPreAction(Action action) {
        actionSet.preAction = action;
    }
    
    /**
     * Sets an action to be performed after the primary action executed
     * by the device.  This is useful for modifying the behavior of concrete
     * subclasses of the device without accessing their primary action.
     */
    public void setPostAction(Action action) {
        actionSet.postAction = action;
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
    public Dimension getDimension() {return unit.dimension();} //may want to override this in most cases
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
    
    private static class ActionSet implements Action, java.io.Serializable {
        Action preAction, postAction, action;
        
        public void actionPerformed() {
            if(preAction != null) preAction.actionPerformed();
            action.actionPerformed();
            if(postAction != null) postAction.actionPerformed();
        }
        
        public String getLabel() {
            return action.getLabel();
        }
    }
}
