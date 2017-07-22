/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Fraction;

import java.io.ObjectStreamException;

/**
 * Decimal representation of something that represents the fractional 
 * amount of a whole (e.g., mole fraction) as a decimal value typically
 * between 0 and 1.
 */
public class Decimal extends SimpleUnit {

  /**
   * Singleton instance of this unit 
   */
	public static final Decimal UNIT = new Decimal();

	private Decimal() {
       super(Fraction.DIMENSION,
        	1.0,
        	"Decimal",
        	"",
        	Prefix.NOT_ALLOWED
        	);
	}
}
