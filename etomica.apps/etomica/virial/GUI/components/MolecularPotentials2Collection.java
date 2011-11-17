package etomica.virial.GUI.components;

import java.util.HashMap;

import etomica.api.IAtomType;
import etomica.api.IPotential;
import etomica.potential.PotentialGroup;

public class MolecularPotentials2Collection extends PotentialCollections{

	
	//Array Object to declare Molecular Potential - Pairwise Interactions
		private IPotential Species1Molecular;
		private IPotential Species2Molecular;
		
		private PotentialGroup pInterGroupIJ;
		
		public PotentialGroup getpInterGroupIJ() {
			return pInterGroupIJ;
		}
		public void setpInterGroupIJ(PotentialGroup pInterGroupIJ) {
			this.pInterGroupIJ = pInterGroupIJ;
		}

		//HashMap to store unlike-Atom Sets
		private	HashMap<String[],IAtomType[]> AtomSetMix;
		
		//HashMap to store all Potentials Sets common for like & unlike interactions
		private	HashMap<String[],IPotential> PotentialSets;
		
		
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
		public IPotential getSpecies1Molecular() {
			return Species1Molecular;
		}
		public void setSpecies1Molecular(IPotential species1Molecular) {
			Species1Molecular = species1Molecular;
		}
		public IPotential getSpecies2Molecular() {
			return Species2Molecular;
		}
		public void setSpecies2Molecular(IPotential species2Molecular) {
			Species2Molecular = species2Molecular;
		}
		
		
		//Methods from Abstract Class
		
		public void setAtomSetsMix(HashMap<String[],IAtomType[]> AtomSet){
			setAtomSetMix(AtomSet);
		}
		
		public void setMolecularPotentialPure(int i, IPotential MPotential){
			switch(i){
			case 1: setSpecies1Molecular(MPotential);
					break;
			case 2: setSpecies2Molecular(MPotential);
					break;
			default:setSpecies1Molecular(MPotential);
					break;		
			}

		}
		
		public void setPotentialSets(HashMap<String[],IPotential> PotentialSets){
			setPotentialSetsAtomic(PotentialSets);
		}
		
		public void setInterPotentialGroupIJ(PotentialGroup pGroup){
			setpInterGroupIJ(pGroup);
		}
		
		public HashMap<String[],IAtomType[]> getAtomSetsMix(){
			return getAtomSetMix();}
		
		public IPotential getMolecularPotentialPure(int i){
			IPotential TempPotential;
			switch(i){
			case 1: TempPotential = getSpecies1Molecular();
					break;
			case 2: TempPotential = getSpecies2Molecular();
					break;
			default:TempPotential = getSpecies1Molecular();
					break;		
			}
			return TempPotential;

		}
		
		public HashMap<String[],IPotential> getPotentialSets(){
			return getPotentialSetsAtomic();}
		
		public PotentialGroup getInterPotentialGroupIJ(){
			return getpInterGroupIJ();
		}
		
		
}
