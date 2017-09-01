/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import etomica.units.Prefix;
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

import java.io.ObjectStreamException;

/**
 * Dimension for all power units. Simulation unit of power is D-A^2/ps^3
 */
public final class Power extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Power();
    /**
     * The simulation unit of power is D-A^2/ps^3.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1.0, "sim power units", "D-A^2/ps^3", Prefix.NOT_ALLOWED);

    private Power() {
        super("Power", 2, 1, -3);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.power();
    }

}
