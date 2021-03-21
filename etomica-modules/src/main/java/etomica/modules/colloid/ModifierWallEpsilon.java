/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.colloid;

import etomica.modifier.Modifier;
import etomica.potential.P1HardFieldGeneric;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;

/**
 * Acts as a modifier for the square-well epsilon for use with P2HardGeneric
 */
public class ModifierWallEpsilon implements Modifier {

    private final P1HardFieldGeneric p;

    public ModifierWallEpsilon(P1HardFieldGeneric p) {
        this.p = p;
    }

    @Override
    public void setValue(double newValue) {
        p.setEnergy(0, -newValue);
        p.setEnergy(2, -newValue);
    }

    @Override
    public double getValue() {
        return -p.getEnergyForState(0);
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
