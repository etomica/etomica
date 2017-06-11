/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Frequency;

/**
 * The Hertz unit of frequency, equal to 1/second.
 */
public final class Hertz extends SimpleUnit {

  /**
   * Singleton instance of this class.
   */
    public static final Hertz UNIT = new Hertz();

    public Hertz() {
        super(Frequency.DIMENSION,
        	1e-12, //conversion from 1/s to 1/ps
        	"Hertz", "Hz", Prefix.ALLOWED
        	);   
    }

}
