/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.ActionGroupSeries;
import etomica.action.ActionToggle;
import etomica.action.IAction;
import etomica.action.controller.Controller;
import etomica.modifier.ModifierBoolean;

/**
 * Button that toggles a boolean value. This device can connect to any object
 * capable of switching between two states. The device operates through a
 * ModifierBoolean instance that must be connected to the state of the
 * controlled object. Action is performed through the controller if given at
 * construction or set afterwards.
 * 
 * @author David Kofke
 */
public class DeviceToggleButton extends DeviceButton {

    public DeviceToggleButton(Controller controller) {
        super(controller);
    }
    
    public DeviceToggleButton(Controller controller, ModifierBoolean modifier) {
        this(controller, modifier, "True", "False");
    }

    public DeviceToggleButton(Controller controller, ModifierBoolean modifier,
                              String trueText, String falseText) {
        super(controller);
        setModifier(modifier, trueText, falseText);
    }

    public void setModifier(ModifierBoolean modifier, String trueText, String falseText) {
        toggleAction = new ActionToggle(modifier, trueText, falseText);
        setAction(new ActionGroupSeries(toggleAction, relabelButton));
        button.setText(modifier.getBoolean() ? trueText : falseText);
    }

    /**
     * Toggles the button to the given state. No action is performed if the
     * button is already in the given state.
     */
    public void setState(boolean b) {
        if (b != getState())
            button.doClick();
    }

    /**
     * Returns the current true/false state of the button.
     */
    public boolean getState() {
        return toggleAction.getModifier().getBoolean();
    }

    /**
     * Specifies the button's label when the toggle is set to true.
     */
    public void setTrueLabel(String text) {
        toggleAction.setTrueLabel(text);
        if (getState())
            button.setText(getTrueLabel());
    }

    /**
     * @return the label displayed when the toggle is set to true.
     */
    public String getTrueLabel() {
        return toggleAction.getTrueLabel();
    }

    /**
     * Specifies the button's label when the toggle is set to false.
     */
    public void setFalseLabel(String text) {
        toggleAction.setFalseLabel(text);
        if (!getState())
            button.setText(getFalseLabel());
    }

    /**
     * @return the label displayed when the toggle is set to false.
     */
    public String getFalseLabel() {
        return toggleAction.getFalseLabel();
    }

    protected ActionToggle toggleAction;
    private final IAction relabelButton = new IAction() {

        public void actionPerformed() {
            DeviceToggleButton.this.setLabel(toggleAction.getLabel());
        }
    };
}
