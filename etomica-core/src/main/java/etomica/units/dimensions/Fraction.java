/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import java.io.ObjectStreamException;

import etomica.units.Decimal;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Dimension for a quantity representing a fractional amount.  Examples
 * include a decimal, percent, parts-per-million, etc.
 */
public final class Fraction extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = Null.DIMENSION;
    /**
     * The simulation unit is the Decimal.
     */
    public static final Unit SIM_UNIT = Decimal.UNIT;

    private Fraction() {
        super("Fraction", 0, 0, 0);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.fraction();
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
