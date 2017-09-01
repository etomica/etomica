/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modifier;

import etomica.units.dimensions.Dimension;

/**
 * A Modifier object permits changes to be made to a property of another object.
 * The target property or object may be unknown to the object that is using the modifier.
 * For example, a modifier may be constructed to adjust a certain property, and then
 * the modifier can be passed to a generic Device, thereby connecting the device to the property.
 * @author David Kofke
 */
public interface Modifier {
    
    /**
     * Sets the value of the property.
     */
    public abstract void setValue(double newValue);
    
    /**
     * Gets the current value of the property.
     * @return
     */
    public abstract double getValue();
    
    /**
     * Returns the dimension of the property being modified.
     */
    public abstract Dimension getDimension();

    /**
     * Returns a descriptive label of the property being modified.
     */
    public abstract String getLabel();
    
}
