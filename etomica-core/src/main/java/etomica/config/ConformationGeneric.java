/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtomList;
import etomica.space.Vector;

public class ConformationGeneric implements IConformation {

	protected final Vector[] coords;
	
	public ConformationGeneric(Vector[] coords) {
		this.coords = coords;
	}
	
	public void initializePositions(IAtomList atomList) {
		for (int i = 0; i<atomList.size(); i++) {
			atomList.get(i).getPosition().E(coords[i]);
		}
	}

}
