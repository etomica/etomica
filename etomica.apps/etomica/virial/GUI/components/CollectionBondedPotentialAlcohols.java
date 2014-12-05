/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;


import java.util.HashMap;

import etomica.api.IPotential;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;

import etomica.potential.P3BondAngle;

import etomica.space.ISpace;
import etomica.units.Kelvin;



public class CollectionBondedPotentialAlcohols extends ACollectionBondedPotential{
	
	
	
	@Override
	public void addBondedPotentialSets(CollectionPotentialAtomicLike pObject,
			ISpace space, int speciesIndex) {
		// TODO Auto-generated method stub
		
		pObject.setHashMapPotentialIntraBonded(new HashMap<String[],IPotential>());
		pObject.setHashMapAtomsetIterator(new HashMap<String[],AtomsetIteratorBasisDependent>());
		
		P3BondAngle P3 = new P3BondAngle(space);
        P3.setAngle(108.5*Math.PI/180);
        P3.setEpsilon(Kelvin.UNIT.toSim(55400));
        int[][] triplets = new int[][] {{0,1,2}};
        pObject.getHashMapPotentialIntraBonded().put(new String[]{"P3","P3"}, P3);
        pObject.getHashmapAtomsetIterator().put(new String[]{"P3","P3"}, new Atomset3IteratorIndexList(triplets));
	}

}
