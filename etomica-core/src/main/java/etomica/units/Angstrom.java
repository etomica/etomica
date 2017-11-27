/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Length;

import java.io.ObjectStreamException;

public final class Angstrom extends SimpleUnit {

    /**
     * Singleton instance of this unit
     */
    public static final Angstrom UNIT = new Angstrom();

    private Angstrom() {
        super(Length.DIMENSION, 1.0,// conversion to simulation units
                "angstroms", "\u00c5", // unicode for the Angstrom symbol
                Prefix.NOT_ALLOWED);
    }
}
