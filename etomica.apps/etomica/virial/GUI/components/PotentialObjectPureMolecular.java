package etomica.virial.GUI.components;

import etomica.api.IPotential;

public class PotentialObjectPureMolecular extends PotentialObject{
	
	//Array Object to declare Molecular Potential - Pairwise Interactions
	private IPotential Species1Molecular;

	public IPotential getSpecies1Molecular() {
		return Species1Molecular;
	}

	public void setSpecies1Molecular(IPotential species1Molecular) {
		Species1Molecular = species1Molecular;
	}
	
	//Methods from abstarct class
	public void setMolecularPotentialPure(IPotential MPotential){
		setSpecies1Molecular(MPotential);

	}
	
	public IPotential getMolecularPotentialPure(){
		return getSpecies1Molecular();
	}
	
	
	

}
