/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.views;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.Border;

import etomica.virial.GUI.models.ModelSimulationEnvironment;

public class ViewSimEnvironmentSelection extends JFrame{
	
	private JLabel message;
	
	private JButton closeWindow;
	private JButton saveValues;
	
	private JTextField temperatureField;
	private JTextField noOfStepsField;
	
	private JTextField[] sigmaHSRefField;
	private JTextField sigmaHSRefFieldA;
	private JTextField sigmaHSRefFieldB;
	private JTextField sigmaHSRefFieldC;
	
	private JLabel[] sigmaHSRefLabel;
	private Border[] sigmaHSRefBorder;
	private JPanel SimEnvParamPanel;
	
	public ViewSimEnvironmentSelection(ModelSimulationEnvironment simEnvironmentObj){
			super();
			this.setMinimumSize(new Dimension(500,400));
			
			SimEnvParamPanel = new JPanel();
			GridBagLayout gridbagOtherParam = new GridBagLayout();
			SimEnvParamPanel.setLayout(gridbagOtherParam);
		    this.add(SimEnvParamPanel);
			
			closeWindow = new JButton("Cancel");
			saveValues = new JButton("Save");
			
		    temperatureField = new JTextField();
			temperatureField.setColumns(10);
			temperatureField.setSize(new Dimension((int)temperatureField.getSize().getWidth(),20));
			temperatureField.setText(Double.toString(simEnvironmentObj.getTemperature()));
			JLabel temperatureLabel = new JLabel("Temperature(K)");
			temperatureLabel.setLabelFor(temperatureField);
			
			Border compoundField1;
			compoundField1 = BorderFactory.createCompoundBorder(
		    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
			temperatureField.setBorder(compoundField1);

			
			noOfStepsField = new JTextField();
			noOfStepsField.setColumns(10);
			JLabel NoOfStepsLabel = new JLabel("Steps");
			noOfStepsField.setText(Integer.toString(simEnvironmentObj.getNoOfSteps()));
			NoOfStepsLabel.setLabelFor(noOfStepsField);
			
			
			Border compoundField2;
			compoundField2 = BorderFactory.createCompoundBorder(
		    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
			noOfStepsField.setBorder(compoundField2);
			
			sigmaHSRefField = new JTextField[10];
			sigmaHSRefLabel = new JLabel[10];
			sigmaHSRefBorder = new Border[10];
			for(int i=0;i<10;i++){
				sigmaHSRefField[i] = new JTextField();
				sigmaHSRefField[i].setColumns(10);
				String index = Integer.toString(i+1);
				sigmaHSRefLabel[i] = new JLabel("SigmaHSRef"+index);
				sigmaHSRefLabel[i].setLabelFor(sigmaHSRefField[i]);
				if(simEnvironmentObj.getSigmaHSRef(i) == 0.0){
					sigmaHSRefField[i].setText("0.0");
					sigmaHSRefField[i].setVisible(false);
					sigmaHSRefLabel[i].setVisible(false);
				}else{
					sigmaHSRefField[i].setText(Double.toString(simEnvironmentObj.getSigmaHSRef(i)));
				}
				
				sigmaHSRefBorder[i] = BorderFactory.createCompoundBorder(
			    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
				sigmaHSRefField[i].setBorder(sigmaHSRefBorder[i]);
			}
			

			
			JComponent[] OtherParamLeft= {temperatureLabel, NoOfStepsLabel, sigmaHSRefLabel[0],sigmaHSRefLabel[1],sigmaHSRefLabel[2],sigmaHSRefLabel[3],
											sigmaHSRefLabel[4],sigmaHSRefLabel[5],sigmaHSRefLabel[6],sigmaHSRefLabel[7],
												sigmaHSRefLabel[8],sigmaHSRefLabel[9],saveValues};
			JComponent[] OtherParamRight = { temperatureField,noOfStepsField,sigmaHSRefField[0],sigmaHSRefField[1],sigmaHSRefField[2],sigmaHSRefField[3],
												sigmaHSRefField[4],sigmaHSRefField[5],sigmaHSRefField[6],sigmaHSRefField[7],
													sigmaHSRefField[8],sigmaHSRefField[8],closeWindow};
		    addLeftRightComponents(OtherParamLeft,OtherParamRight,gridbagOtherParam,SimEnvParamPanel);
		    

		    Border compound2;
		    compound2 = BorderFactory.createCompoundBorder(
		    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
		    SimEnvParamPanel.setBorder(compound2);

	}
	
	public JTextField getTemperatureField() {
		return temperatureField;
	}

	public JTextField getNoOfStepsField() {
		return noOfStepsField;
	}

	public JTextField getSigmaHSRefField(int index) {
		return sigmaHSRefField[index];
	}
	
	public void setSigmaHSRefField(JTextField textField, int index) {
		sigmaHSRefField[index] = textField;
	}
	
	public JTextField getSigmaHSRefFieldA() {
		return sigmaHSRefFieldA;
	}
	
	public JTextField getSigmaHSRefFieldB() {
		return sigmaHSRefFieldB;
	}
	
	public JTextField getSigmaHSRefFieldC() {
		return sigmaHSRefFieldC;
	}

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
	
	
	public void addChooseOkayButtonListener(ActionListener Comp) {
		closeWindow.addActionListener(Comp);
	}

	
	public JButton getCloseWindow() {
		return closeWindow;
	}

	public JButton getSaveValues() {
		return saveValues;
	}
	
	public void openSimEnvParamView(){
		this.setVisible(true);
	}
	
	public void closeSimEnvParamView(){
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

}
