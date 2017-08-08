/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modifier;

import etomica.action.IAction;


/**
 * Wraps a ModifierBoolean instance so that the modification can be performed
 * as the consequence of an action. Permits for example a modification
 * initiated by a user-interface thread to be executed by the controller
 * thread.  Call to setValueForAction stores value for modifier, which
 * is then invoked when a subsequent call to actionPerformed is made.
 *
 * @author David Kofke, Andrew Schultz
 */
public class ModifyBooleanAction implements IAction {

    /**
     * Constuctor requires the wrapped ModifierBoolean
     */
    public ModifyBooleanAction(ModifierBoolean modifier) {
        this.modifier = modifier;
    }

    public void actionPerformed() {
        modifier.setBoolean(value);
    }

    /**
     * Stores the new value to be set by the wrapped modifier.  Value is
     * not set until actionPerformed is invoked.
     */
    public void setValueForAction(boolean newBoolean) {
        value = newBoolean;
    }

    public ModifierBoolean getWrappedModifier() {
        return modifier;
    }

    private final ModifierBoolean modifier;
    private boolean value;
}
