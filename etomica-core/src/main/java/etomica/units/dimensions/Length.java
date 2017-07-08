/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import java.io.ObjectStreamException;

import etomica.units.Angstrom;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Dimension for all units of length. Simulation unit of length is the angstrom.
 */
public final class Length extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Length();
    /**
     * Simulation unit for length is the angstrom.
     */
    public static final Unit SIM_UNIT = Angstrom.UNIT;
 
    private Length() {
        super("Length", 1, 0, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.length();
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
