/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modifier;

import etomica.potential.P2HardGeneric;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;

/**
 * Acts as a modifier for the square-well epsilon for use with P2HardGeneric
 */
public class ModifierSQWEpsilon implements Modifier {

    private final P2HardGeneric p;

    public ModifierSQWEpsilon(P2HardGeneric p) {
        this.p = p;
    }

    @Override
    public void setValue(double newValue) {
        p.setEnergyForState(1, -newValue);
    }

    @Override
    public double getValue() {
        return -p.getEnergyForState(1);
    }

    @Override
    public Dimension getDimension() {
        return Energy.DIMENSION;
    }

    @Override
    public String getLabel() {
        return "SQW epsilon";
    }
}
