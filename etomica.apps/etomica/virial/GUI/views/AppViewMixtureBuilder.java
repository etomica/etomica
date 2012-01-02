package etomica.virial.GUI.views;

import etomica.virial.GUI.containers.AppFrameMixtureBuilder;
import etomica.virial.GUI.models.AppModelMixtureBuilder;


public class AppViewMixtureBuilder {
		
	
	private ViewSpeciesSelection viewSpeciesSelection;
	
	

	public AppViewMixtureBuilder(AppFrameMixtureBuilder mixtureBuilder, AppModelMixtureBuilder SM) {
		
		viewSpeciesSelection = new ViewSpeciesSelection(SM.getModelSpeciesSelection(),SM.getModelSelectedSpecies());
	    mixtureBuilder.add(viewSpeciesSelection);
	    mixtureBuilder.setFrameProperties();
	}
	
	public ViewSpeciesSelection getViewSpeciesSelection() {
		return viewSpeciesSelection;
	}

	public void setViewSpeciesSelection(ViewSpeciesSelection viewSpeciesSelection) {
		this.viewSpeciesSelection = viewSpeciesSelection;
	}
}
