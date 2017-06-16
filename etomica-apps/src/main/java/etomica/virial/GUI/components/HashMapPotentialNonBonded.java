/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;

import etomica.potential.IPotential;

import java.util.HashMap;

public class HashMapPotentialNonBonded {

	public HashMap<String[],IPotential> hashMapPotentialNonBonded;
	
	
	private HashMapPotentialNonBonded(){
		hashMapPotentialNonBonded = new HashMap<String[],IPotential>();
	}
	
	
	private static class HolderHashMapPotentialNonBonded { 
         public static final HashMapPotentialNonBonded instance = new HashMapPotentialNonBonded();
	}

	public static HashMapPotentialNonBonded getInstance() {
         return HolderHashMapPotentialNonBonded.instance;
	}

	public HashMap<String[], IPotential> getHashMapPotentialNonBonded() {
		return hashMapPotentialNonBonded;
	}

	public void setHashMapPotentialNonBonded(
			HashMap<String[], IPotential> hashMapPotentialNonBonded) {
		this.hashMapPotentialNonBonded = hashMapPotentialNonBonded;
	}
	
	
	
	
}
