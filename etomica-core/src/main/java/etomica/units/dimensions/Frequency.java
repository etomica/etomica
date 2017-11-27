/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import etomica.units.*;
import etomica.units.systems.UnitSystem;

/**
 * Dimension for all units of time frequency, 1/Time
 */
public final class Frequency extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Frequency();
    /**
     * The simulation unit is 1/ps = 1/(10^-12 s) = 10^12/s = 1 THz.
     */
    public static final Unit SIM_UNIT = new PrefixedUnit(Prefix.TERA,Hertz.UNIT);

    private Frequency() {
        super("Frequency", 0, 0, -1);// LMTCtNl;
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return new CompoundUnit(new Unit[] {unitSystem.time()},new double[] {0, 0,-1});
    }

}
