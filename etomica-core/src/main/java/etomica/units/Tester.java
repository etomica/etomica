/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Force;

public class Tester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		UnitGraphics gUI = new UnitGraphics();
		Dimension targetDimension = Force.DIMENSION;
		gUI.startWithDim(targetDimension);
	}
}
