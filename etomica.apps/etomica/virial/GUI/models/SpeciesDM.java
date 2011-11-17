package etomica.virial.GUI.models;

public class SpeciesDM {
	
	private CreateSpeciesDM_IFactory species1DataModel;
	private CreateSpeciesDM_IFactory species2DataModel;
	
	private int species1MoleculeCount;
	private int species2MoleculeCount;
	
	//AddFlag in Species1
	private boolean species1Added;
	
	//ChangesSpecies in Species1
	private boolean speciesSelectionChanged;
	
	public SpeciesDM(){
		reset();
	}
	
	public CreateSpeciesDM_IFactory getSpecies1DataModel() {
		return species1DataModel;
	}

	public void setSpecies1DataModel(CreateSpeciesDM_IFactory species1DataModel) {
		this.species1DataModel = species1DataModel;
	}

	public CreateSpeciesDM_IFactory getSpecies2DataModel() {
		return species2DataModel;
	}

	public void setSpecies2DataModel(CreateSpeciesDM_IFactory species2DataModel) {
		this.species2DataModel = species2DataModel;
	}

	public int getSpecies1MoleculeCount() {
		return species1MoleculeCount;
	}

	public void setSpecies1MoleculeCount(int species1MoleculeCount) {
		this.species1MoleculeCount = species1MoleculeCount;
	}

	public int getSpecies2MoleculeCount() {
		return species2MoleculeCount;
	}

	public void setSpecies2MoleculeCount(int species2MoleculeCount) {
		this.species2MoleculeCount = species2MoleculeCount;
	}

	public boolean isSpecies1Added() {
		return species1Added;
	}

	public void setSpecies1Added(boolean species1Added) {
		this.species1Added = species1Added;
	}

	public boolean isSpeciesSelectionChanged() {
		return speciesSelectionChanged;
	}

	public void setSpeciesSelectionChanged(boolean speciesSelectionChanged) {
		this.speciesSelectionChanged = speciesSelectionChanged;
	}

	/** Reset to initial value. */
    public void reset() {
    	species1DataModel = null;
    	species1DataModel = null;
    	species1MoleculeCount = 0;
    	species2MoleculeCount = 0;
    	speciesSelectionChanged = false;
    	species1Added = false;
    }
}
