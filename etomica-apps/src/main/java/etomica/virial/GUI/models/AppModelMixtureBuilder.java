/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
