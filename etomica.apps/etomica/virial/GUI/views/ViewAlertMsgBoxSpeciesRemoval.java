/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.views;

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class ViewAlertMsgBoxSpeciesRemoval extends JFrame{
	
	/**
	 * 
	 */
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JLabel Message;
	
	private JButton yesRemoveSpeciesButton;
	private JButton noDoNotRemoveSpeciesButton;

	public ViewAlertMsgBoxSpeciesRemoval(String message){
		//JFrame MessageFrame = new JFrame("Alert");
		super();
		this.setMinimumSize(new Dimension(400,100));
		this.setBounds(700, 200, 0, 0);
		
		Message = new JLabel(message);
		yesRemoveSpeciesButton = new JButton("Yes");
		noDoNotRemoveSpeciesButton = new JButton("No");
		
		this.setLayout(new FlowLayout());
		this.add(Message);
		this.add(yesRemoveSpeciesButton);
		this.add(noDoNotRemoveSpeciesButton);
		this.setVisible(true);

	}

	public JButton getYesButton() {
		return yesRemoveSpeciesButton;
	}

	public void setYesButton(JButton closeWindow) {
		yesRemoveSpeciesButton = closeWindow;
	}
	
	public void addYesRemoveSpeciesButtonListener(ActionListener mal){
		yesRemoveSpeciesButton.addActionListener(mal);
	}
	
	public void addNoDoNotRemoveSpeciesButtonListener(ActionListener mal){
		noDoNotRemoveSpeciesButton.addActionListener(mal);
	}
}
