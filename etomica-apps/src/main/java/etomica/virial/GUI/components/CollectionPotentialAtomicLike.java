/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import etomica.atom.IAtomType;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.potential.PotentialGroup;
import etomica.virial.MCMoveClusterTorsionMulti;

public class CollectionPotentialAtomicLike implements ICollectionPotential{

	private int speciesIndex;
	
	private PotentialGroup potentialGroupInterNonBondedLike;
	
	private PotentialGroup potentialGroupIntraBonded;

	private MCMoveClusterTorsionMulti[] mcMoveClusterTorsionMulti;
	
	//Array of HashMap to store both like-Molecule AtomSets
	private HashMap<String[],IAtomType[]> hashMapAtomTypesLikePairs;

	
	//HashMap to store all Potentials Sets common for like & unlike interactions
	private	HashMapPotentialNonBonded hashmapPotentialNonBonded;
	
	//HashMap to store all Bonded potentials for species 1 and species 2
	private HashMap<String[],IPotential> hashmapPotentialIntraBonded;
	
	private HashMap<String[],AtomsetIteratorBasisDependent> hashMapAtomsetIterators;
	
	private IPotentialAtomic potentialIntraBondedTorsion;
	
	public CollectionPotentialAtomicLike(int index){
		this.speciesIndex = index;
		setHashMapPotentialNonBonded();
	}
	
	public void setHashMapAtomTypeLikePairs(HashMap<String[], IAtomType[]> hashMapAtomTypesLikePairs) {
		// TODO Auto-generated method stub
		this.hashMapAtomTypesLikePairs = hashMapAtomTypesLikePairs;
		
	}

	
	public void setHashMapPotentialNonBonded() {
		// TODO Auto-generated method stub
		this.hashmapPotentialNonBonded = HashMapPotentialNonBonded.getInstance();
	}

	
	public void setHashMapPotentialIntraBonded(
			HashMap<String[], IPotential> hashMapPotentialIntraBonded) {
		this.hashmapPotentialIntraBonded = hashMapPotentialIntraBonded;
		
	}

	
	public void setMCMoveClusterTorsionMulti(
			MCMoveClusterTorsionMulti[] mcMoveClusterTorsionMulti) {
		// TODO Auto-generated method stub
		this.mcMoveClusterTorsionMulti = mcMoveClusterTorsionMulti;
	}

	
	public void setHashMapAtomsetIterator(
			HashMap<String[], AtomsetIteratorBasisDependent> hashMapAtomsetIterators) {
		// TODO Auto-generated method stub
		this.hashMapAtomsetIterators = hashMapAtomsetIterators;
	}

	
	public void setPotentialGroupInterNonBondedLike(PotentialGroup potentialGroupInterNonBondedLike) {
		// TODO Auto-generated method stub
		this.potentialGroupInterNonBondedLike = potentialGroupInterNonBondedLike;
	}

	
	public void setPotentialGroupIntraBonded(PotentialGroup potentialGroupIntraBonded) {
		// TODO Auto-generated method stub
		this.potentialGroupIntraBonded = potentialGroupIntraBonded;
	}

	
	public void setPotentialIntraBondedTorsion(
			IPotentialAtomic potentialIntraBondedTorsion) {
		// TODO Auto-generated method stub
		this.potentialIntraBondedTorsion = potentialIntraBondedTorsion;
	}


	
	public HashMapPotentialNonBonded getHashMapPotentialNonBonded() {
		// TODO Auto-generated method stub
		return this.hashmapPotentialNonBonded;
	}

	
	public HashMap<String[], IAtomType[]> getHashMapAtomTypeLikePairs() {
		// TODO Auto-generated method stub
		return this.hashMapAtomTypesLikePairs;
	}

	
	public HashMap<String[], IPotential> getHashMapPotentialIntraBonded() {
		// TODO Auto-generated method stub
		return this.hashmapPotentialIntraBonded;
	}

	
	public MCMoveClusterTorsionMulti[] getMCMoveClusterTorsionMulti() {
		// TODO Auto-generated method stub
		return this.mcMoveClusterTorsionMulti;
	}

	
	public HashMap<String[], AtomsetIteratorBasisDependent> getHashmapAtomsetIterator() {
		// TODO Auto-generated method stub
		return this.hashMapAtomsetIterators;
	}

	
	public PotentialGroup getPotentialGroupInterNonBondedLike() {
		// TODO Auto-generated method stub
		return this.potentialGroupInterNonBondedLike;
	}

	
	public PotentialGroup getPotentialGroupIntraBonded() {
		// TODO Auto-generated method stub
		return this.potentialGroupIntraBonded;
	}

	
	public IPotentialAtomic getPotentialIntraBondedTorsion() {
		// TODO Auto-generated method stub
		return this.potentialIntraBondedTorsion;
	}

	public void addToHashMapPotentialNonBonded(String[] hashKey, IPotential hashValuePotential){
		Map<String[], IPotential> nonBondedPotentialsMap = this.hashmapPotentialNonBonded.getHashMapPotentialNonBonded();
		Set<Entry<String[], IPotential>> nonBondedPotentialEntries = nonBondedPotentialsMap.entrySet();
		Iterator<Entry<String[], IPotential>> nonBondedPotentialsItr = nonBondedPotentialEntries.iterator();
		
		if(nonBondedPotentialsMap.size() == 0){
			this.hashmapPotentialNonBonded.getHashMapPotentialNonBonded().put(hashKey,hashValuePotential);
		}else{
			int tableIterationIndex = 0;
			while(nonBondedPotentialsItr.hasNext()){
				Map.Entry nonBondedPotentialsEntry = nonBondedPotentialsItr.next();
				String[] nonBondedPotentialsMapKey= (String[]) nonBondedPotentialsEntry.getKey();
				if(!(nonBondedPotentialsMapKey[0] == hashKey[0] && nonBondedPotentialsMapKey[1] == hashKey[1]) &&
						!(nonBondedPotentialsMapKey[0] == hashKey[1] && nonBondedPotentialsMapKey[1] == hashKey[0])){
					
					tableIterationIndex++;
				}
			}
			if(tableIterationIndex == nonBondedPotentialsMap.size()){
				this.hashmapPotentialNonBonded.getHashMapPotentialNonBonded().put(hashKey, hashValuePotential);
			}
			
		}
	}

	public int getSpeciesIndex() {
		return speciesIndex;
	}

	public void setSpeciesIndex(int speciesIndex) {
		this.speciesIndex = speciesIndex;
	}
	
	
	
}
