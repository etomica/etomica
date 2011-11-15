package etomica.virial.GUI.components;

import java.util.HashMap;

import etomica.api.IAtomType;
import etomica.api.IPotential;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.potential.PotentialGroup;
import etomica.virial.MCMoveClusterTorsionMulti;

public class AtomicPotentialCollection extends PotentialCollectionFactory{

	private PotentialGroup pInterGroupII;
	
	private PotentialGroup pIntraGroupII;

	

	private MCMoveClusterTorsionMulti[] torsionMovesI;
	
	//Array of HashMap to store both like-Molecule AtomSets
	private HashMap<String[],IAtomType[]> AtomSetPure;

	
	//HashMap to store all Potentials Sets common for like & unlike interactions
	private	HashMap<String[],IPotential> PotentialSets;
	
	//HashMap to store all Bonded potentials for species 1 and species 2
	private HashMap<String[],IPotential> BondedIIPotentialSets;
	
	private HashMap<String[],AtomsetIteratorBasisDependent> IteratorSets;
	
	public PotentialGroup getpIntraGroupII() {
		return pIntraGroupII;
	}

	public void setpIntraGroupII(PotentialGroup pIntraGroupII) {
		this.pIntraGroupII = pIntraGroupII;
	}

	public PotentialGroup getpInterGroupII() {
		return pInterGroupII;
	}

	public void setpInterGroupII(PotentialGroup pInterGroupII) {
		this.pInterGroupII = pInterGroupII;
	}

	public MCMoveClusterTorsionMulti[] getTorsionMovesI() {
		return torsionMovesI;
	}

	public void setTorsionMovesI(MCMoveClusterTorsionMulti[] torsionMovesI) {
		this.torsionMovesI = torsionMovesI;
	}

	public HashMap<String[], IAtomType[]> getAtomSetPure() {
		return AtomSetPure;
	}

	public void setAtomSetPure(HashMap<String[], IAtomType[]> atomSetPure) {
		AtomSetPure = atomSetPure;
	}

	public HashMap<String[], IPotential> getPotentialSetsAtomic() {
		return PotentialSets;
	}

	public void setPotentialSetsAtomic(HashMap<String[], IPotential> potentialSets) {
		PotentialSets = potentialSets;
	}

	public HashMap<String[], IPotential> getBondedIIPotentialSets() {
		return BondedIIPotentialSets;
	}

	public void setBondedIIPotentialSets(
			HashMap<String[], IPotential> bondedIIPotentialSets) {
		BondedIIPotentialSets = bondedIIPotentialSets;
	}

	public HashMap<String[], AtomsetIteratorBasisDependent> getIteratorSetsAtomic() {
		return IteratorSets;
	}

	public void setIteratorSetsAtomic(
			HashMap<String[], AtomsetIteratorBasisDependent> iteratorSets) {
		IteratorSets = iteratorSets;
	}
	
	//Methods from Abstract Class
	
	public void setInterPotentialGroupII(PotentialGroup pGroup){
		setpInterGroupII(pGroup);
	}

	public void setIntraPotentialGroup(PotentialGroup pGroup){
		setpIntraGroupII(pGroup);
				
	}
	

	public void setAtomSetsPure(HashMap<String[],IAtomType[]> AtomSet){
		setAtomSetPure(AtomSet);
	}
	
	public void setPotentialSets(HashMap<String[],IPotential> PotentialSets){
		setPotentialSetsAtomic(PotentialSets);
	}
	

	public void setBondedPotentialSets(HashMap<String[],IPotential> BondedPotentialSets){
		setBondedIIPotentialSets(BondedPotentialSets);
	}
	
	
	public void setMCMoveClusterTorsionMulti(MCMoveClusterTorsionMulti[] TorsionMoves){
		setTorsionMovesI(TorsionMoves);
	}
	
	
	public void setIteratorSets(HashMap<String[],AtomsetIteratorBasisDependent> IteratorSets){
		setIteratorSetsAtomic(IteratorSets);
	}
	
	
	public PotentialGroup getInterPotentialGroupII(){
		return getpInterGroupII();}
	
	public PotentialGroup getIntraPotentialGroup(){
		return getpIntraGroupII();}
	
	public HashMap<String[],IAtomType[]> getAtomSetsPure(){
		return getAtomSetPure();}
	
	public HashMap<String[],IPotential> getPotentialSets(){
		return getPotentialSetsAtomic();}
	
	public HashMap<String[],IPotential> getBondedPotentialSets(){
		return getBondedIIPotentialSets();
	}
	
	public MCMoveClusterTorsionMulti[] getMCMoveClusterTorsionMulti(){
		return getTorsionMovesI();
	};
	
	public HashMap<String[],AtomsetIteratorBasisDependent> getIteratorSets(){
		return getIteratorSetsAtomic();
	};
	
	
		
}
