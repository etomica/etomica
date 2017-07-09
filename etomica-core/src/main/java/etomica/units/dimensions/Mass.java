/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import java.io.ObjectStreamException;

import etomica.units.Dalton;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Dimension for all units of mass. Simulation unit for mass is the Dalton (or atomic
 * mass unit, 1/N_Avo grams).
 */
public final class Mass extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Mass();
    /**
     * The simulation unit is the Dalton.
     */
    public static final Unit SIM_UNIT = Dalton.UNIT;

    private Mass() {
        super("Mass", 0, 1, 0);// LMTCtNl;
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.mass();
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
