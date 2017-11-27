/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import java.io.ObjectStreamException;

import etomica.units.Picosecond;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Dimension for all units of time.
 */
public final class Time extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Time();
    /**
     * Simulation unit for time is the picosecond.
     */
    public static final Unit SIM_UNIT = Picosecond.UNIT;

    private Time() {
        super("Time", 0, 0, 1);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.time();
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
