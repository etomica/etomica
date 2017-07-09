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
 * Undefined dimension used for quantities with undefined or unknown dimensions.
 */
public final class Undefined extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Undefined();
    /**
     * Unit for a quantity having undefined or unknown units. 
     * Any conversion using this unit will cause output to be NaN (not a number).
     */
    public static final Unit UNIT = new SimpleUnit(DIMENSION, Double.NaN, "undefined", "", Prefix.NOT_ALLOWED);

    /**
     * Private constructor for singleton instantiation.
     */
    private Undefined() {
        super("Undefined",
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
    }
    
    /**
     * Returns true if the given object is the singleton instance of this class.
     */
    public boolean equals(Object obj) {
        return obj == this;
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
