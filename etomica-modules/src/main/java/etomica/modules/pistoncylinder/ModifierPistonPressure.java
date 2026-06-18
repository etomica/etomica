/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.modifier.Modifier;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;

/**
 * Modifier for the piston pressure.  Should be followed by a IntegratorPistonUpdate action.
 */
public class ModifierPistonPressure implements Modifier {
    public ModifierPistonPressure(Space _space, P1HardMovingBoundary potential,
                                  Dimension pressureDimension) {
        pistonPotential = potential;
        pressureDim = pressureDimension;
        space = _space;
    }

    public void setValue(double p) {
        pistonPotential.setPressure(space.D() == 3 ? -p : p);
    }

    public double getValue() {
        double p = pistonPotential.getPressure();
        return space.D() == 3 ? -p : p;
    }

    public Dimension getDimension() {
        return pressureDim;
    }

    public String getLabel() {
        return "Piston pressure";
    }

    public String toString() {
        return getLabel();
    }

    private final P1HardMovingBoundary pistonPotential;
    private final Dimension pressureDim;
    private Space space;
}
