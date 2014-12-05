/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.controllers;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import etomica.virial.GUI.models.ModelAlertMsgDialog;
import etomica.virial.GUI.models.ModelSpeciesSelection;
import etomica.virial.GUI.views.ViewAlertMsgBox;
import etomica.virial.GUI.views.ViewAlertMsgBoxSpeciesRemoval;
import etomica.virial.GUI.views.ViewSpeciesSelection;

public class ControllerAlertMsgBoxSpeciesRemoval {

	private int speciesIndex;
	private ModelAlertMsgDialog alertMsgBoxDM;
	private ModelSpeciesSelection modelSpeciesSelection;


	private ViewAlertMsgBoxSpeciesRemoval alertMsgBoxViewSpeciesRemoval;

	
	public ControllerAlertMsgBoxSpeciesRemoval(String AlertMessage, ModelSpeciesSelection modelSpeciesSelection,int index){
		this.modelSpeciesSelection = modelSpeciesSelection;

		speciesIndex = index;
		this.alertMsgBoxDM = new ModelAlertMsgDialog(AlertMessage);
		alertMsgBoxViewSpeciesRemoval = new ViewAlertMsgBoxSpeciesRemoval(alertMsgBoxDM.getAlertMessage());
		alertMsgBoxViewSpeciesRemoval.addYesRemoveSpeciesButtonListener(new YesRemoveSpeciesButtonListener());
		alertMsgBoxViewSpeciesRemoval.addNoDoNotRemoveSpeciesButtonListener(new NoDoNotRemoveSpeciesButtonListener());
		
		
	}
	

	public class YesRemoveSpeciesButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			alertMsgBoxViewSpeciesRemoval.setVisible(false);
			modelSpeciesSelection.setRemoveSpeciesAfterAlertFlag(true);
			//Change code here!!! to Reset the potential2 choice!!!!
			//modelSpeciesSelection.reset();
		}
	 
	 }
	
	public class NoDoNotRemoveSpeciesButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			alertMsgBoxViewSpeciesRemoval.setVisible(false);
			modelSpeciesSelection.setRemoveSpeciesAfterAlertFlag(false);
		}
	 
	 }
	
}

