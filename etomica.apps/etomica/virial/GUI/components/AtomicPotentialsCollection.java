package etomica.virial.GUI.components;

import java.util.HashMap;

import etomica.api.IAtomType;
import etomica.api.IPotential;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.potential.PotentialGroup;
import etomica.virial.MCMoveClusterTorsionMulti;

public class AtomicPotentialsCollection extends PotentialCollections {
	
	
	
	
	
	
	private PotentialGroup pInterGroupII;
	private PotentialGroup pInterGroupIJ;
	private PotentialGroup pInterGroupJJ;
	
	//Array Object to declare all others...
	private PotentialGroup pIntraGroupII;
	private PotentialGroup pIntraGroupJJ;
	
	private MCMoveClusterTorsionMulti[] torsionMovesI;
	private MCMoveClusterTorsionMulti[] torsionMovesJ;
	
	//Array of HashMap to store both like-Molecule AtomSets
	private HashMap<String[],IAtomType[]> AtomSetPureII;
	private HashMap<String[],IAtomType[]> AtomSetPureJJ;
						
	//HashMap to store unlike-Atom Sets
	private	HashMap<String[],IAtomType[]> AtomSetMix;
						
	//HashMap to store all Potentials Sets common for like & unlike interactions
	private	HashMap<String[],IPotential> PotentialSets;
		
		
	//HashMap to store all Bonded potentials for species 1 and species 2
	private HashMap<String[],IPotential> BondedIIPotentialSetsII;
	private HashMap<String[],IPotential> BondedIIPotentialSetsJJ;
		

	private HashMap<String[],AtomsetIteratorBasisDependent> IteratorSetsII;
	private HashMap<String[],AtomsetIteratorBasisDependent> IteratorSetsJJ;


	//Getter Setter methods for specific class
	public PotentialGroup getpInterGroupII() {
		return pInterGroupII;
	}


	public void setpInterGroupII(PotentialGroup pInterGroupII) {
		this.pInterGroupII = pInterGroupII;
	}


	public PotentialGroup getpInterGroupIJ() {
		return pInterGroupIJ;
	}


	public void setpInterGroupIJ(PotentialGroup pInterGroupIJ) {
		this.pInterGroupIJ = pInterGroupIJ;
	}


	public PotentialGroup getpInterGroupJJ() {
		return pInterGroupJJ;
	}


	public void setpInterGroupJJ(PotentialGroup pInterGroupJJ) {
		this.pInterGroupJJ = pInterGroupJJ;
	}


	public PotentialGroup getpIntraGroupII() {
		return pIntraGroupII;
	}


	public void setpIntraGroupII(PotentialGroup pIntraGroupII) {
		this.pIntraGroupII = pIntraGroupII;
	}


	public PotentialGroup getpIntraGroupJJ() {
		return pIntraGroupJJ;
	}


	public void setpIntraGroupJJ(PotentialGroup pIntraGroupJJ) {
		this.pIntraGroupJJ = pIntraGroupJJ;
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


	public HashMap<String[], IAtomType[]> getAtomSetPureII() {
		return AtomSetPureII;
	}


	public void setAtomSetPureII(HashMap<String[], IAtomType[]> atomSetPureII) {
		AtomSetPureII = atomSetPureII;
	}


	public HashMap<String[], IAtomType[]> getAtomSetPureJJ() {
		return AtomSetPureJJ;
	}


	public void setAtomSetPureJJ(HashMap<String[], IAtomType[]> atomSetPureJJ) {
		AtomSetPureJJ = atomSetPureJJ;
	}


	public HashMap<String[], IAtomType[]> getAtomSetMix() {
		return AtomSetMix;
	}


	public void setAtomSetMix(HashMap<String[], IAtomType[]> atomSetMix) {
		AtomSetMix = atomSetMix;
	}


	public HashMap<String[], IPotential> getPotentialSetsAtomic() {
		return PotentialSets;
	}


	public void setPotentialSetsAtomic(HashMap<String[], IPotential> potentialSets) {
		PotentialSets = potentialSets;
	}


	public HashMap<String[], IPotential> getBondedIIPotentialSetsII() {
		return BondedIIPotentialSetsII;
	}


	public void setBondedIIPotentialSetsII(
			HashMap<String[], IPotential> bondedIIPotentialSetsII) {
		BondedIIPotentialSetsII = bondedIIPotentialSetsII;
	}


	public HashMap<String[], IPotential> getBondedIIPotentialSetsJJ() {
		return BondedIIPotentialSetsJJ;
	}


	public void setBondedIIPotentialSetsJJ(
			HashMap<String[], IPotential> bondedIIPotentialSetsJJ) {
		BondedIIPotentialSetsJJ = bondedIIPotentialSetsJJ;
	}


	public HashMap<String[], AtomsetIteratorBasisDependent> getIteratorSetsII() {
		return IteratorSetsII;
	}


	public void setIteratorSetsII(
			HashMap<String[], AtomsetIteratorBasisDependent> iteratorSetsII) {
		IteratorSetsII = iteratorSetsII;
	}


	public HashMap<String[], AtomsetIteratorBasisDependent> getIteratorSetsJJ() {
		return IteratorSetsJJ;
	}


	public void setIteratorSetsJJ(
			HashMap<String[], AtomsetIteratorBasisDependent> iteratorSetsJJ) {
		IteratorSetsJJ = iteratorSetsJJ;
	}


	//Methods from Abstract Class
	
	public void setAtomSetsPure(int i, HashMap<String[],IAtomType[]> AtomSet){
		switch(i){
		case 1: setAtomSetPureII(AtomSet);
				break;
		case 2: setAtomSetPureJJ(AtomSet);
				break;
		default:setAtomSetPureII(AtomSet);
				break;
		}
	}
	
	public void setAtomSetsMix(HashMap<String[],IAtomType[]> AtomSet){
		setAtomSetMix(AtomSet);
	}
	
	public void setPotentialSets(HashMap<String[],IPotential> PotentialSets){
		setPotentialSetsAtomic(PotentialSets);
	}
	
	
	public void setBondedPotentialSets(int i,HashMap<String[],IPotential> BondedPotentialSets){
		switch(i){
		case 1: setBondedIIPotentialSetsII(BondedPotentialSets);
				break;
		case 2: setBondedIIPotentialSetsJJ(BondedPotentialSets);
				break;
		default:setBondedIIPotentialSetsII(BondedPotentialSets);
				break;
		}
	}
	
	
	public void setMCMoveClusterTorsionMulti(int i,MCMoveClusterTorsionMulti[] TorsionMoves){
		switch(i){
		case 1: setTorsionMovesI(TorsionMoves);
				break;
		case 2: setTorsionMovesJ(TorsionMoves);
				break;
		default:setTorsionMovesI(TorsionMoves);
				break;
		}
	}
	
	
	public void setIteratorSets(int i,HashMap<String[],AtomsetIteratorBasisDependent> IteratorSets){
		switch(i){
		case 1: setIteratorSetsII(IteratorSets);
				break;
		case 2: setIteratorSetsJJ(IteratorSets);
				break;
		default:setIteratorSetsII(IteratorSets);
				break;
		}
	}
	
	public void setInterPotentialGroupII(int i, PotentialGroup pGroup){
		switch(i){
		case 1: setpInterGroupII(pGroup);
				break;
		case 2: setpInterGroupJJ(pGroup);
				break;
		default:setpInterGroupII(pGroup);
				break;		
		}
	}
	

	public void setInterPotentialGroupIJ(PotentialGroup pGroup){
		setpInterGroupIJ(pGroup);
	}
	
	public void setIntraPotentialGroup(int i, PotentialGroup pGroup){
		switch(i){
		case 1: setpIntraGroupII(pGroup);
				break;
		case 2: setpIntraGroupJJ(pGroup);
				break;
		default:setpIntraGroupII(pGroup);
				break;		
		}
	}
	
	
	public HashMap<String[],IAtomType[]> getAtomSetsPure(int i){
		HashMap<String[],IAtomType[]> TempMap;
		switch(i){
		case 1: TempMap = getAtomSetPureII();
				break;
		case 2: TempMap = getAtomSetPureJJ();
				break;
		default:TempMap = getAtomSetPureII();
				break;
		}
		return TempMap;
	}
	
	public HashMap<String[],IAtomType[]> getAtomSetsMix(){
		return getAtomSetMix();}
	
	public HashMap<String[],IPotential> getPotentialSets(){
		return getPotentialSetsAtomic();}
	
	public HashMap<String[], IPotential> getBondedPotentialSets(int i){
		HashMap<String[], IPotential> TempPotential;
		switch(i){
		case 1: TempPotential = getBondedIIPotentialSetsII();
				break;
		case 2: TempPotential = getBondedIIPotentialSetsJJ();
				break;
		default:TempPotential = getBondedIIPotentialSetsII();
				break;
		}
		return TempPotential;
	}
	
	
	public MCMoveClusterTorsionMulti[] getMCMoveClusterTorsionMulti(int i){
		MCMoveClusterTorsionMulti[] TempTorsionMoves;
		switch(i){
		case 1: TempTorsionMoves = getTorsionMovesI();
				break;
		case 2: TempTorsionMoves = getTorsionMovesJ();
				break;
		default:TempTorsionMoves = getTorsionMovesI();
				break;
		}
		return TempTorsionMoves;
	}
	

	public HashMap<String[],AtomsetIteratorBasisDependent> getIteratorSets(int i){
		HashMap<String[],AtomsetIteratorBasisDependent> TempIterator;
		switch(i){
		case 1: TempIterator = getIteratorSetsII();
				break;
		case 2: TempIterator = getIteratorSetsJJ();
				break;
		default:TempIterator = getIteratorSetsII();
				break;
		}
		return TempIterator;
	}
	
	public PotentialGroup getInterPotentialGroupII(int i){
		PotentialGroup TempGroup;
		switch(i){
		case 1: TempGroup = getpInterGroupII();
				break;
		case 2: TempGroup = getpInterGroupJJ();
				break;
		default:TempGroup = getpInterGroupII();
				break;
		}
		return TempGroup;
	}
	
	public PotentialGroup getInterPotentialGroupIJ(){
		return getpInterGroupIJ();
	}
	
	public PotentialGroup getIntraPotentialGroupII(int i){
		PotentialGroup TempGroup;
		switch(i){
		case 1: TempGroup = getpIntraGroupII();
				break;
		case 2: TempGroup = getpIntraGroupJJ();
				break;
		default:TempGroup = getpIntraGroupII();
				break;
		}
		return TempGroup;
	}
	
	
}
