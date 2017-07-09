/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import java.io.ObjectStreamException;

import etomica.units.Prefix;
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Dimension specified dimensionless quantities which have no other interpretation
 * (e.g., the quantity is not known to be an angle, or a fraction).
 */
public final class Null extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Null();
    
    /**
     * Unit for an otherwise unspecified dimensionless quantity.
     */
    public static final Unit UNIT = new SimpleUnit(DIMENSION, 1.0, "", "", Prefix.NOT_ALLOWED);

    /**
     * Private constructor for singleton instantiation.
     */
    private Null() {
        super("Dimensionless", 0, 0, 0);
    }

    /**
     * Returns UNIT, regardless of the given unit system.
     */
    public Unit getUnit(UnitSystem unitSystem) {
        return UNIT;
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
