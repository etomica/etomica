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

public class ControllerAlertMsgBox {

	private ModelAlertMsgDialog alertMsgBoxDM;
	private ModelSpeciesSelection modelSpeciesSelection;

	private ViewAlertMsgBox alertMsgBoxView;

	
	public ControllerAlertMsgBox(String AlertMessage, ModelSpeciesSelection modelSpeciesSelection){
		this.modelSpeciesSelection = modelSpeciesSelection;
		this.alertMsgBoxDM = new ModelAlertMsgDialog(AlertMessage);
		
			alertMsgBoxView = new ViewAlertMsgBox(alertMsgBoxDM.getAlertMessage());
			alertMsgBoxView.addCloseWindowButtonListener(new CloseWindowButtonListener());
		
		
	}
	
	public class CloseWindowButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			alertMsgBoxView.setVisible(false);
			
			//Change code here!!! to Reset the potential2 choice!!!!
			//modelSpeciesSelection.reset();
		}
	 
	 }
	
	
}

