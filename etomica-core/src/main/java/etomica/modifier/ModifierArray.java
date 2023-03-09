/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modifier;

import etomica.units.dimensions.Dimension;

/**
 * Wraps an array of modifiers when use of ModifierGeneral is not possible.
 */
public class ModifierArray implements Modifier {
    private final Modifier[] modifiers;

    public ModifierArray(Modifier[] modifiers) {
        this.modifiers = modifiers;
    }

    @Override
    public void setValue(double newValue) {
        for (Modifier modifier : modifiers) {
            modifier.setValue(newValue);
        }
    }

    @Override
    public double getValue() {
        return modifiers[0].getValue();
    }

    @Override
    public Dimension getDimension() {
        return modifiers[0].getDimension();
    }

    @Override
    public String getLabel() {
        return modifiers[0].getLabel();
    }
}
