package etomica.virial.GUI.components;

import java.util.HashMap;

import etomica.api.IAtomType;
import etomica.api.IPotential;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.potential.PotentialGroup;
import etomica.virial.MCMoveClusterTorsionMulti;

public class MixedPotentialsCollection extends PotentialCollectionFactory{

	
	private PotentialGroup pInterGroupII;
	private PotentialGroup pInterGroupIJ;
	
	private PotentialGroup pIntraGroupII;

	private MCMoveClusterTorsionMulti[] torsionMovesI;
	
	//Array of HashMap to store both like-Molecule AtomSets
	private HashMap<String[],IAtomType[]> AtomSetPureII;
	
	//HashMap to store unlike-Atom Sets
	private	HashMap<String[],IAtomType[]> AtomSetMix;
	
	//HashMap to store all Potentials Sets common for like & unlike interactions
	private	HashMap<String[],IPotential> PotentialSets;
	
	//HashMap to store all Bonded potentials for species 1 and species 2
	private HashMap<String[],IPotential> BondedIIPotentialSetsII;
	
	private HashMap<String[],AtomsetIteratorBasisDependent> IteratorSets;
	
	private IPotential SpeciesMolecular;

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

	public PotentialGroup getpIntraGroupII() {
		return pIntraGroupII;
	}

	public void setpIntraGroupII(PotentialGroup pIntraGroupII) {
		this.pIntraGroupII = pIntraGroupII;
	}

	public MCMoveClusterTorsionMulti[] getTorsionMovesI() {
		return torsionMovesI;
	}

	public void setTorsionMovesI(MCMoveClusterTorsionMulti[] torsionMovesI) {
		this.torsionMovesI = torsionMovesI;
	}

	public HashMap<String[], IAtomType[]> getAtomSetPureII() {
		return AtomSetPureII;
	}

	public void setAtomSetPureII(HashMap<String[], IAtomType[]> atomSetPureII) {
		AtomSetPureII = atomSetPureII;
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

	public HashMap<String[], AtomsetIteratorBasisDependent> getIteratorSetsII() {
		return IteratorSets;
	}

	public void setIteratorSetsII(
			HashMap<String[], AtomsetIteratorBasisDependent> iteratorSets) {
		IteratorSets = iteratorSets;
	}

	public IPotential getSpeciesMolecular() {
		return SpeciesMolecular;
	}

	public void setSpeciesMolecular(IPotential speciesMolecular) {
		SpeciesMolecular = speciesMolecular;
	}

	//Methods from Abstract Class
	
		public void setAtomSetsPure(HashMap<String[],IAtomType[]> AtomSet){
			setAtomSetPureII(AtomSet);
		}
		
		public void setAtomSetsMix(HashMap<String[],IAtomType[]> AtomSet){
			setAtomSetMix(AtomSet);
		}
		
		public void setPotentialSets(HashMap<String[],IPotential> PotentialSets){
			setPotentialSetsAtomic(PotentialSets);
		}
		
		
		public void setBondedPotentialSets(HashMap<String[],IPotential> BondedPotentialSets){
			setBondedIIPotentialSetsII(BondedPotentialSets);
		}
		
		
		public void setMCMoveClusterTorsionMulti(MCMoveClusterTorsionMulti[] TorsionMoves){
			setTorsionMovesI(TorsionMoves);
		}
		
		
		public void setIteratorSets(int i,HashMap<String[],AtomsetIteratorBasisDependent> IteratorSets){
			setIteratorSetsII(IteratorSets);
		}
		
		public void setInterPotentialGroupII(PotentialGroup pGroup){
			setpInterGroupII(pGroup);
					
		}
		public void setInterPotentialGroupIJ(PotentialGroup pGroup){
			setpInterGroupIJ(pGroup);
		}

		public void setIntraPotentialGroup(PotentialGroup pGroup){
			setpIntraGroupII(pGroup);
					
		}
		
		public void setMolecularPotentialPure(IPotential MPotential){
			setSpeciesMolecular(MPotential);

		}
		
		
		//Getter Methods
		
		public HashMap<String[],IAtomType[]> getAtomSetsPure(){
			return getAtomSetPureII();}
		
		
		public HashMap<String[],IAtomType[]> getAtomSetsMix(){
			return getAtomSetMix();}
		

		public HashMap<String[],IPotential> getPotentialSets(){
			return getPotentialSetsAtomic();}
		
		public HashMap<String[],IPotential> getBondedPotentialSets(){
			return getBondedIIPotentialSetsII();
		}
		
		public MCMoveClusterTorsionMulti[] getMCMoveClusterTorsionMulti(){
			return getTorsionMovesI();
		};
		
		public HashMap<String[],AtomsetIteratorBasisDependent> getIteratorSets(){
			return getIteratorSetsII();
		};
		
		public PotentialGroup getInterPotentialGroupII(){
			return getpInterGroupII();
		}
		
		public PotentialGroup getInterPotentialGroupIJ(){
			return getpInterGroupIJ();
		}
		
		public PotentialGroup getIntraPotentialGroup(){
			return getpIntraGroupII();
		}
		
		public IPotential getMolecularPotentialPure(){
			return getSpeciesMolecular();
		}
		
		
		
		
		
		
		
}
