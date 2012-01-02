package etomica.virial.GUI.controllers;

import etomica.virial.GUI.models.AppModelMixtureBuilder;
import etomica.virial.GUI.views.AppViewMixtureBuilder;

public class AppControllerMixtureBuilder {
	
	public ControllerSpeciesSelection controllerSpeciesSelection;
	
	public AppControllerMixtureBuilder(
			AppViewMixtureBuilder appViewMixtureBuilder,
			AppModelMixtureBuilder appModelMixtureBuilder) {
		
		controllerSpeciesSelection = new ControllerSpeciesSelection(appViewMixtureBuilder.getViewSpeciesSelection(), 
																		appModelMixtureBuilder.getModelSpeciesSelection(),
																			appModelMixtureBuilder.getModelSelectedSpecies());
		// TODO Auto-generated constructor stub
	}

}
