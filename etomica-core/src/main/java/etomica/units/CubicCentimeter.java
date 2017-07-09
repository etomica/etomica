/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Volume;

import java.io.ObjectStreamException;

/**
 * The cubic centimeter unit of volume, cm^3.
 */
public final class CubicCentimeter extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final CubicCentimeter UNIT = new CubicCentimeter();

    private CubicCentimeter() {
        super(Volume.DIMENSION, 1e+24, // conversion from cm^3 to Angstroms^3
                "cubic centimeters", "cc", Prefix.NOT_ALLOWED);
    }

    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton UNIT
     */
    private Object readResolve() throws ObjectStreamException {
        return UNIT;
    }

    private static final long serialVersionUID = 1;

}
