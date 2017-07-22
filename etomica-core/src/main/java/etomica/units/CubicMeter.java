/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Volume;

import java.io.ObjectStreamException;

/**
 * The cubic meter unit of volume, m^3.
 */
public final class CubicMeter extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final CubicMeter UNIT = new CubicMeter();

    private CubicMeter() {
        super(Volume.DIMENSION, 1e+30, // conversion from meters^3 to Angstroms^3
                "cubic meters", "m^3", Prefix.NOT_ALLOWED);
    }
}
