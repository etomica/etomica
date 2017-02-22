/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modifier;

/**
 * Interface that permits changes to be made to a boolean property of another object.
 * This type of modifier might be plugged into a checkbox or toggle button.
 */

public interface ModifierBoolean {
    
    /**
     * Sets the modified value to the given boolean.
     */
    public void setBoolean(boolean b);

    /**
     * Returns the current state of the modified value.
     */
    public boolean getBoolean();

}