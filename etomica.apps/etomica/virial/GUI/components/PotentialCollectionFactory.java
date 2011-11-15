package etomica.virial.GUI.components;

import java.util.HashMap;

import etomica.api.IAtomType;
import etomica.api.IPotential;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.potential.PotentialGroup;
import etomica.virial.MCMoveClusterTorsionMulti;

public abstract class PotentialCollectionFactory {
	
	//Specific to Atomic Potentials
	
	public void setAtomSetsPure(int i, HashMap<String[],IAtomType[]> AtomSet){};
	
	public void setAtomSetsPure(HashMap<String[],IAtomType[]> AtomSet){};
	
	
	public void setAtomSetsMix(HashMap<String[],IAtomType[]> AtomSet){};
	
	
	public void setPotentialSets(HashMap<String[],IPotential> PotentialSets){};
	
	
	public void setBondedPotentialSets(int i,HashMap<String[],IPotential> BondedPotentialSets){};
	
	public void setBondedPotentialSets(HashMap<String[],IPotential> BondedPotentialSets){};
	
	
	public void setMCMoveClusterTorsionMulti(int i,MCMoveClusterTorsionMulti[] TorsionMoves){};
	
	public void setMCMoveClusterTorsionMulti(MCMoveClusterTorsionMulti[] TorsionMoves){};
	
	
	public void setIteratorSets(int i,HashMap<String[],AtomsetIteratorBasisDependent> IteratorSets){};
	
	public void setIteratorSets(HashMap<String[],AtomsetIteratorBasisDependent> IteratorSets){};
	
	
	public void setInterPotentialGroupII(int i, PotentialGroup pGroup){};
	
	public void setInterPotentialGroupII(PotentialGroup pGroup){};
	
	public void setInterPotentialGroupIJ(PotentialGroup pGroup){};
	
	
	public void setIntraPotentialGroup(int i, PotentialGroup pGroup){};
	
	public void setIntraPotentialGroup(PotentialGroup pGroup){};
	
	
	//Specific to Molecular Potentials
	
	public void setMolecularPotentialPure(int i, IPotential MPotential){};
	
	public void setMolecularPotentialPure(IPotential MPotential){};
	
	public void setMolecularPotentialCross(IPotential MPotential){};
	
	
	
	
	//Getter Methods
	
	public HashMap<String[],IAtomType[]> getAtomSetsPure(int i){
		return null;};
	
	public HashMap<String[],IAtomType[]> getAtomSetsPure(){
		return null;};
		
	public HashMap<String[],IAtomType[]> getAtomSetsMix(){
		return null;};
	
	public HashMap<String[],IPotential> getPotentialSets(){
		return null;};
		
		
		
	public HashMap<String[],IPotential> getBondedPotentialSets(int i){
		return null;
	};
		
	public HashMap<String[],IPotential> getBondedPotentialSets(){
		return null;
	};
	
	public MCMoveClusterTorsionMulti[] getMCMoveClusterTorsionMulti(int i){
		return null;
	};
	
	public MCMoveClusterTorsionMulti[] getMCMoveClusterTorsionMulti(){
		return null;
	};
	
	public HashMap<String[],AtomsetIteratorBasisDependent> getIteratorSets(int i){
		return null;
	};
	
	public HashMap<String[],AtomsetIteratorBasisDependent> getIteratorSets(){
		return null;
	};
	
	
	public PotentialGroup getInterPotentialGroupII(int i){
		return null;};
	
	public PotentialGroup getInterPotentialGroupII(){
		return null;};
	
	public PotentialGroup getInterPotentialGroupIJ(){
		return null;};
		
	
	public PotentialGroup getIntraPotentialGroup(int i){
		return null;};
	
	public PotentialGroup getIntraPotentialGroup(){
		return null;};
	
		
	
	public IPotential getMolecularPotentialPure(int i){
		return null;};
		
	public IPotential getMolecularPotentialPure(){
		return null;};	
		
	public IPotential getMolecularPotentialCross(){
		return null;};
		
	
	

}
