package etomica.virial.GUI.components;

import etomica.api.IPotential;

public class MolecularPotentialsCollection extends PotentialCollections{

	
	//Array Object to declare Molecular Potential - Pairwise Interactions
		private IPotential Species1Molecular;
		private IPotential Species2Molecular;
		private IPotential CrossPotential;
		
		
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
		public IPotential getCrossPotential() {
			return CrossPotential;
		}
		public void setCrossPotential(IPotential crossPotential) {
			CrossPotential = crossPotential;
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
		
		public void setMolecularPotentialCross(IPotential MPotential){
			setCrossPotential(MPotential);
		}
		
		public IPotential getMolecularPotentialPure(int i){
			IPotential tempPotential;
			switch(i){
			case 1: tempPotential = getSpecies1Molecular();
					break;
			case 2: tempPotential = getSpecies2Molecular();
					break;
			default:tempPotential = getSpecies1Molecular();
					break;		
			}
			return tempPotential;
		}
		
		public IPotential getMolecularPotentialCross(){
			return getCrossPotential();
		}
		
}
