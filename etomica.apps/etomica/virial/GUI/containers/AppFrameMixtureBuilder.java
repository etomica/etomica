/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.containers;

import java.awt.BorderLayout;
import java.awt.Dimension;

import javax.swing.JFrame;



public class AppFrameMixtureBuilder extends JFrame{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	
	public AppFrameMixtureBuilder(String FrameLabel){
		super(FrameLabel);
		
	}
	
	

	public void setFrameProperties(){
		this.setVisible(true);
		this.setResizable(true);
		
		//this.setExtendedState(java.awt.Frame.MAXIMIZED_BOTH);  
		this.pack();
	}


	
}
