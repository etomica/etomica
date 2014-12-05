/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

public class ModelAlertMsgDialog {
	
	private String alertMessage;
	
	private ModelSpeciesSelection modelSpeciesSelection;
	
	public ModelAlertMsgDialog(String alertMessage){
		setAlertMessage(alertMessage);
	}

	public ModelSpeciesSelection getSpeciesSelectionDM() {
		return modelSpeciesSelection;
	}

	public void setSpeciesSelectionDM(ModelSpeciesSelection modelSpeciesSelection) {
		this.modelSpeciesSelection = modelSpeciesSelection;
	}

	public String getAlertMessage() {
		return this.alertMessage;
	}

	public void setAlertMessage(String alertMessage) {
		this.alertMessage = alertMessage;
	}
	
}
