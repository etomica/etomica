package etomica.virial.GUI.components;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import etomica.api.IPotential;

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
