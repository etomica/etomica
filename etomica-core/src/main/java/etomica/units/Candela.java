/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.LuminousIntensity;

import java.io.ObjectStreamException;

/**
 * The candela is the luminous intensity, in a given direction, of a source that
 * emits monochromatic radiation of frequency 540 x 1012 hertz and that has a
 * radiant intensity in that direction of 1/683 watt per steradian.
 */
public final class Candela extends SimpleUnit {

    /**
     * Singleton instance of this unit.
     */
    public static final Candela UNIT = new Candela();

    private Candela() {
        super(LuminousIntensity.DIMENSION, 1.0, "Candela", "cd", Prefix.ALLOWED);
    }

}
