/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Quantity;

import java.io.ObjectStreamException;

/**
 * Unit class for the number, or quantity, of something. This is one of the
 * default simulation units.
 */
public final class Count extends SimpleUnit {

    /**
     * Singleton instance of this unit
     */
    public static final Count UNIT = new Count();

    private Count() {
        super(Quantity.DIMENSION, 1.0, "Count", "", Prefix.NOT_ALLOWED);
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
