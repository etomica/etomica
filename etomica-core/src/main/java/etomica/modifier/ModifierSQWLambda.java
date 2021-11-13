/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modifier;

import etomica.potential.P2HardGeneric;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Null;

/**
 * Acts as a modifier for the square-well lambda for use with P2HardGeneric
 */
public class ModifierSQWLambda implements Modifier {

    private final P2HardGeneric p;

    public ModifierSQWLambda(P2HardGeneric p) {
        this.p = p;
    }

    @Override
    public void setValue(double newValue) {
        p.setCollisionDiameter(1, newValue * p.getCollisionDiameter(0));
    }

    @Override
    public double getValue() {
        return p.getCollisionDiameter(1) / p.getCollisionDiameter(0);
    }

    @Override
    public Dimension getDimension() {
        return Null.DIMENSION;
    }

    @Override
    public String getLabel() {
        return "SQW lambda";
    }
}
