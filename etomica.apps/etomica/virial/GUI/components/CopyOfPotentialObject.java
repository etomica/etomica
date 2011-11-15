package etomica.virial.GUI.components;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import etomica.api.IAtomType;
import etomica.api.IPotential;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.potential.PotentialGroup;
import etomica.virial.MCMoveClusterTorsionMulti;

public class CopyOfPotentialObject {

	//Array Object to into Account For Allowing like-like paiwise interactions- Atomic Potential Interactions
	//For eg. long chain alkanes, we ll include any intercations beyond the  
	private PotentialGroup[] pInterGroupII;
	private PotentialGroup pInterGroupIJ;
	
	//Array Object to declare all others...
	private PotentialGroup[] pIntraGroupII;
	
	private MCMoveClusterTorsionMulti[] torsionMovesI;
	private MCMoveClusterTorsionMulti[] torsionMovesJ;
	
	

	//Array Object to declare Molecular Potential - Pairwise Interactions
	private IPotential Species1Molecular;
	private IPotential Species2Molecular;
	private IPotential CrossPotential;
	

	//Array of HashMap to store both like-Molecule AtomSets
	private HashMap<String[],IAtomType[]>[] AtomSetPure;
					
	//HashMap to store unlike-Atom Sets
	private	HashMap<String[],IAtomType[]> AtomSetMix;
					
	//HashMap to store all Potentials Sets common for like & unlike interactions
	private	HashMap<String[],IPotential> PotentialSets;
	
	
	//HashMap to store all Bonded potentials for species 1 and species 2
	private HashMap<String[],IPotential>[] BondedIIPotentialSets;
	

	private HashMap<String[],AtomsetIteratorBasisDependent>[] IteratorSets;
	
	
	
	

	CopyOfPotentialObject(){
		pInterGroupII = new PotentialGroup[2];
		pInterGroupIJ = new PotentialGroup(2);
		pIntraGroupII = new PotentialGroup[2];
		
		AtomSetPure = new HashMap[2];
		AtomSetMix = new HashMap<String[],IAtomType[]>();
		PotentialSets = new HashMap<String[],IPotential>();
		BondedIIPotentialSets = new HashMap[2];
		IteratorSets = new HashMap[2];
	}

	
	//pIntergroupII for pure species 1 and pure species 2
	
	public PotentialGroup getpInterGroupII(int i) {
		return pInterGroupII[i];
	}

	public void setpInterGroupII(PotentialGroup pInterGroupII,int i) {
		this.pInterGroupII[i] = pInterGroupII;
	}

	//pIntergroupIJ for potentials formed for unlike pairs
	public PotentialGroup getpInterGroupIJ() {
		return pInterGroupIJ;
	}

	public void setpInterGroupIJ(PotentialGroup pInterGroupIJ) {
		this.pInterGroupIJ = pInterGroupIJ;
	}
	
	//pIntra group for pure species 1 and pure species 2

	public PotentialGroup getpIntraGroupII(int i) {
		return pIntraGroupII[i];
	}

	public void setpIntraGroupII(PotentialGroup pIntraGroupII,int i) {
		this.pIntraGroupII[i] = pIntraGroupII;
	}

	
	//Molecular potential stored for Species1
	
	public Object getSpecies1Molecular() {
		return Species1Molecular;
	}

	public void setSpecies1Molecular(IPotential species1Molecular) {
		Species1Molecular = species1Molecular;
	}
	
	//Molecular potential stored for Species2

	public Object getSpecies2Molecular() {
		return Species2Molecular;
	}

	public void setSpecies2Molecular(IPotential species2Molecular) {
		Species2Molecular = species2Molecular;
	}

	
	////Molecular cross potentials
	
	public Object getCrossPotential() {
		return CrossPotential;
	}

	public void setCrossPotential(IPotential crossPotential) {
		CrossPotential = crossPotential;
	}

	
	//AtomSet for species 1 and species 2
	
	public HashMap<String[],IAtomType[]>[] getAtomSetPure() {
		return AtomSetPure;
	}
	
	public HashMap<String[],IAtomType[]> getAtomSetPureSpeciesI(int i) {
		return AtomSetPure[i];
	}
	
	public void setAtomSetPure(HashMap<String[],IAtomType[]>[] atomSetPure) {
		AtomSetPure = atomSetPure;
	}
	
	public void setAtomSetPureSpeciesI(HashMap<String[],IAtomType[]> atomSetPure, int i) {
		AtomSetPure[i] = atomSetPure;
	}

	
	//Atomset for Mix
	
	public HashMap<String[],IAtomType[]> getAtomSetMix() {
		return AtomSetMix;
	}

	public void setAtomSetMix(HashMap<String[],IAtomType[]> atomSetMix) {
		AtomSetMix = atomSetMix;
	}

	
	//All potentials sets for species 1, species 2 and unlike pairs
	public HashMap<String[],IPotential> getPotentialSets() {
		return PotentialSets;
	}

	public void setPotentialSets(HashMap<String[],IPotential> potentialSets) {
		PotentialSets = potentialSets;
	}

	
	//BondedPotentials for Pure species 1 and 2
	public HashMap<String[],IPotential>[] getBondedIIPotentialSets() {
		return BondedIIPotentialSets;
	}

	public void setBondedIIPotentialSets(HashMap<String[],IPotential>[] bondedIIPotentialSets) {
		BondedIIPotentialSets = bondedIIPotentialSets;
	}
	
	public HashMap<String[],IPotential> getBondedIIPotentialSetsSpeciesI(int index) {
		return BondedIIPotentialSets[index];
	}
	
	public void setBondedIIPotentialSetsSpeciesI(HashMap<String[],IPotential> bondedIIPotentialSets, int index) {
		BondedIIPotentialSets[index] = bondedIIPotentialSets;
	}
	
	
	public HashMap<String[],AtomsetIteratorBasisDependent>[] getIteratorSets() {
		return IteratorSets;
	}
	
	public HashMap<String[],AtomsetIteratorBasisDependent> getIteratorSetsSpeciesI(int index) {
		return IteratorSets[index];
	}

	public void setIteratorSets(HashMap<String[],AtomsetIteratorBasisDependent>[] iteratorSets) {
		IteratorSets = iteratorSets;
	}
	
	public void setIteratorSetsSpeciesI(HashMap<String[],AtomsetIteratorBasisDependent> iteratorSets, int index) {
		IteratorSets[index] = iteratorSets;
	}
	
	public MCMoveClusterTorsionMulti[] getTorsionMovesI() {
		return torsionMovesI;
	}

	public void setTorsionMovesI(MCMoveClusterTorsionMulti[] torsionMovesI) {
		this.torsionMovesI = torsionMovesI;
	}

	public MCMoveClusterTorsionMulti[] getTorsionMovesJ() {
		return torsionMovesJ;
	}

	public void setTorsionMovesJ(MCMoveClusterTorsionMulti[] torsionMovesJ) {
		this.torsionMovesJ = torsionMovesJ;
	}
	
	
}
