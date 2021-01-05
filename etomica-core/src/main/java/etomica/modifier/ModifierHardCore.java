/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modifier;

import etomica.potential.P2HardGeneric;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;

/**
 * Acts as a modifier for the hard core diameter for use with P2HardGeneric
 */
public class ModifierHardCore implements Modifier {

    private final P2HardGeneric p;

    public ModifierHardCore(P2HardGeneric p) {
        this.p = p;
    }

    @Override
    public void setValue(double newValue) {
        p.setCollisionDiameter(0, newValue);
    }

    @Override
    public double getValue() {
        return p.getCollisionDiameter(0);
    }

    @Override
    public Dimension getDimension() {
        return Length.DIMENSION;
    }

    @Override
    public String getLabel() {
        return "core diameter";
    }
}
