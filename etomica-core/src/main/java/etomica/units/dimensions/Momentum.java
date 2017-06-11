/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import etomica.units.CompoundUnit;
import etomica.units.Prefix;
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * Dimension for all units of momentum, Mass-Length/Time
 */
public final class Momentum extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Momentum();
    /**
     * The simulation unit is A/ps.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "D-angstrom/ps", "D-\u00c5/ps", Prefix.NOT_ALLOWED);

    private Momentum() {
        super("Momentum", 1, 1, -1);// LMTCtNl;
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return new CompoundUnit(new Unit[] {unitSystem.mass(),unitSystem.length(),unitSystem.time()},new double[] {+1, +1,-1});
    }

}
