/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
