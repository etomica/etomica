/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.views;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;

import etomica.space1d.Space1D;

public class ViewAlkaneLengthSelection  extends JFrame {
	
	/**
	 * 
	 */
	
	
	private static final long serialVersionUID = 1L;
	public JButton closeWindow;
	public JButton saveValues;
	
	public JTextField NoOfSpheres;
	
	public JPanel alkaneSpherePanel;
	//private static final NAlkaneSpheresParameter INSTANCE = new NAlkaneSpheresParameter();
	
	public ViewAlkaneLengthSelection(String frameTitle){
		super(frameTitle);
		this.setMinimumSize(new Dimension(500,200));
		
		alkaneSpherePanel = new JPanel();
		this.add(alkaneSpherePanel);
		closeWindow = new JButton("Cancel");
		saveValues = new JButton("Save");
		GridBagLayout gridbagOtherParam = new GridBagLayout();
		alkaneSpherePanel.setLayout(gridbagOtherParam);
	    
		NoOfSpheres = new JTextField();
		NoOfSpheres.setColumns(10);
		NoOfSpheres.setSize(new Dimension(5,20));
		NoOfSpheres.setText("4");
		JLabel NLabel = new JLabel("Cn (where n is 4 or greater)");
		NLabel.setLabelFor(NoOfSpheres);
		
		
		Border compoundField1;
		compoundField1 = BorderFactory.createCompoundBorder(
	    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
		NoOfSpheres.setBorder(compoundField1);
		
		JComponent[] OtherParamLeft= {NLabel,saveValues};
		JComponent[] OtherParamRight = { NoOfSpheres,closeWindow};
	    addLeftRightComponents(OtherParamLeft,OtherParamRight,gridbagOtherParam,alkaneSpherePanel);
		
		
	}
	/*
	public static NAlkaneSpheresParameter getInstance() {
		return INSTANCE;
		
	}*/
	private void addLeftRightComponents(JComponent[] ComponentLeft,
			JComponent[] ComponentRight,
            GridBagLayout gridbag,
            Container container) {
		
			GridBagConstraints c = new GridBagConstraints();
			c.anchor = GridBagConstraints.NORTHEAST;
			int numLabels = ComponentLeft.length;

			for (int i = 0; i < numLabels; i++) {
				c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
				c.fill = GridBagConstraints.NONE;      //reset to default
				c.weightx = 0.0;   
				c.insets = new Insets(0,10,0,0);//reset to default
				container.add(ComponentLeft[i], c);

				c.gridwidth = GridBagConstraints.REMAINDER;     //end row
				c.fill = GridBagConstraints.HORIZONTAL;
				c.weightx = 0.0;
				c.insets = new Insets(0,10,3,0);
				container.add(ComponentRight[i], c);
			}
	}
	public JTextField getNoOfSpheres() {
		return NoOfSpheres;
	}

	public void setNoOfSpheres(JTextField noOfSpheres) {
		NoOfSpheres = noOfSpheres;
	}

	public JButton getCloseWindow() {
		return closeWindow;
	}

	public void setCloseWindow(JButton closeWindow) {
		this.closeWindow = closeWindow;
	}

	public JButton getSaveValues() {
		return saveValues;
	}

	public void setSaveValues(JButton saveValues) {
		this.saveValues = saveValues;
	}
	
	public void openAlkaneView(){
		this.setVisible(true);
	}
	
	public void closeAlkaneView(){
		this.setVisible(false);
	}
	
	//Listener method to run the simulation when the Run button is pressed
	public void addCloseWindowButtonListener(ActionListener mal) {
		closeWindow.addActionListener(mal);
	}
	
	//Listener method to run the simulation when the Run button is pressed
	public void addSaveValuesButtonListener(ActionListener mal) {
		saveValues.addActionListener(mal);
	}
	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				
				ViewAlkaneLengthSelection s = new ViewAlkaneLengthSelection("Species1 n-alkanes");
				s.openAlkaneView();
				s.closeAlkaneView();
				
				
			}
		});
    }
	
}
