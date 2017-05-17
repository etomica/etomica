/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;

import etomica.space.Space;

public abstract class ACollectionBondedPotential {
	
	
	public abstract void  addBondedPotentialSets(CollectionPotentialAtomicLike pObject, Space space, int speciesIndex);

	

}
