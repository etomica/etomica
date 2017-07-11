/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import etomica.units.CompoundUnit;
import etomica.units.Prefix;
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

import java.io.ObjectStreamException;

/**
 * Dimension for all units of velocity, Length/Time
 */
public final class Velocity extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Velocity();
    /**
     * The simulation unit is A/ps.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "angstrom/ps", "\u00c5/ps", Prefix.NOT_ALLOWED);

    private Velocity() {
        super("Velocity", 1, 0, -1);// LMTCtNl;
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return new CompoundUnit(new Unit[] {unitSystem.length(),unitSystem.time()},new double[] {+1,-1});
    }

}
