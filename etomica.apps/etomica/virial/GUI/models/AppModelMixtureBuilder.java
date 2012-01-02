package etomica.virial.GUI.models;

public class AppModelMixtureBuilder {

	public ModelSpeciesSelection modelSpeciesSelection;
	public ModelSelectedSpecies modelSelectedSpecies;
	
	public AppModelMixtureBuilder(){
		modelSpeciesSelection = new ModelSpeciesSelection();
		modelSelectedSpecies = new ModelSelectedSpecies();
	}

	public ModelSpeciesSelection getModelSpeciesSelection() {
		return modelSpeciesSelection;
	}

	public void setModelSpeciesSelection(ModelSpeciesSelection modelSpeciesSelection) {
		this.modelSpeciesSelection = modelSpeciesSelection;
	}

	public ModelSelectedSpecies getModelSelectedSpecies() {
		return modelSelectedSpecies;
	}

	public void setModelSelectedSpecies(ModelSelectedSpecies modelSelectedSpecies) {
		this.modelSelectedSpecies = modelSelectedSpecies;
	}
	
	
}
