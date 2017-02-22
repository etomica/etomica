/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI;

import javax.swing.SwingUtilities;

import etomica.virial.GUI.containers.AppFrameMixtureBuilder;
import etomica.virial.GUI.controllers.AppControllerMixtureBuilder;
import etomica.virial.GUI.controllers.ControllerSpeciesSelection;

import etomica.virial.GUI.models.AppModelMixtureBuilder;
import etomica.virial.GUI.views.AppViewMixtureBuilder;




public class ApplicationMixtureBuilder{
	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				AppFrameMixtureBuilder appFrameMixtureBuilder = new AppFrameMixtureBuilder("Mixture Builder");
				//MainFrame DefaultValueFrame = new MainFrame("Default Values");
			
				AppModelMixtureBuilder appModelMixtureBuilder = new AppModelMixtureBuilder();
				//Instantiate the view for MixtureParameters
				
				AppViewMixtureBuilder appViewMixtureBuilder = new AppViewMixtureBuilder(appFrameMixtureBuilder,appModelMixtureBuilder);
				
				AppControllerMixtureBuilder appControllerMixtureBuilder = new AppControllerMixtureBuilder(appViewMixtureBuilder,appModelMixtureBuilder);
				
				
				System.out.println("MixtureBuilder ConsoleOutput");
				
			}
		});
    }

}
