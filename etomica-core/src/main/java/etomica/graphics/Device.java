/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.IAction;
import etomica.action.activity.Controller;
import etomica.units.Unit;
import etomica.units.dimensions.Dimension;

import java.awt.*;

/**
 * Base class for all Devices.  These are objects that permit manipulation of the
 * fields of other objects.  They provide a (usually) graphical interface to permit
 * values of fields to be input and changed.
 */
public abstract class Device {

    protected Unit unit;
    protected Controller controller;
    private final ActionSet actionSet = new ActionSet();

    public Device() {
        this(null);
    }

    public Device(Controller controller) {
        this.controller = controller;
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
    protected void doAction(IAction action) {
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
    public void setPreAction(IAction action) {
        actionSet.preAction = action;
    }
    
    /**
     * Returns the action to be performed before the primary action executed
     * by the device.
     */
    public IAction getPreAction() {
        return actionSet.preAction;
    }

    /**
     * Sets an action to be performed after the primary action executed
     * by the device.  This is useful for modifying the behavior of concrete
     * subclasses of the device without accessing their primary action.
     */
    public void setPostAction(IAction action) {
        actionSet.postAction = action;
    }

    public IAction getPostAction() {
        return actionSet.postAction;
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

    public abstract Component graphic();

    public void setUnit(Unit u) {
        unit = u;
    }

    public Unit getUnit() {
        return unit;
    }

    public Dimension getDimension() {
        return unit.dimension();
    } //may want to override this in most cases

    protected static class ActionSet implements IAction {
        IAction preAction, postAction, action;

        public void actionPerformed() {
            if (preAction != null) preAction.actionPerformed();
            action.actionPerformed();
            if(postAction != null) postAction.actionPerformed();
        }
    }
}
