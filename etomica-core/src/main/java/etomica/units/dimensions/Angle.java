/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import java.io.ObjectStreamException;

import etomica.units.Radian;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Base for all angular units. Simulation unit for angles is radians.
 */
public class Angle extends Dimension {

    public static final Dimension DIMENSION = new Angle();
    public static final Unit SIM_UNIT = Radian.UNIT;

    private Angle() {
        super("Angle", 0, 0, 0);
    }

     public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.angle();
    }

    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton DIMENSION
     */
    private Object readResolve() throws ObjectStreamException {
        return DIMENSION;
    }
    
    private static final long serialVersionUID = 1;

}
