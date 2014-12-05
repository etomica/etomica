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
import javax.swing.JTextArea;

public class ViewAlertMsgBox extends JFrame{
	
	private static final long serialVersionUID = 1L;
	private JTextArea Message;
	
	private JButton CloseWindow;

	public ViewAlertMsgBox(String message){
		//JFrame MessageFrame = new JFrame("Alert");
		super();
		this.setMinimumSize(new Dimension(400,100));
		this.setBounds(700, 200, 0, 0);
		
		Message = new JTextArea(message);
		Message.setEditable(false);
		CloseWindow = new JButton("Ok");
		
		this.setLayout(new FlowLayout());
		this.add(Message);
		this.add(CloseWindow);
		this.setVisible(true);

	}

	public JButton getCloseWindow() {
		return CloseWindow;
	}

	public void setCloseWindow(JButton closeWindow) {
		CloseWindow = closeWindow;
	}
	
	public void addCloseWindowButtonListener(ActionListener mal){
		CloseWindow.addActionListener(mal);
	}

	

}
