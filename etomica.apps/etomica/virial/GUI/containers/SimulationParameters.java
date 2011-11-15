package etomica.virial.GUI.containers;

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
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.Border;

import etomica.virial.GUI.components.SimulationEnvironmentObject;

public class SimulationParameters extends JPanel{
	
	private JLabel Message;
	
	private JButton CloseWindow;
	private JButton SaveValues;
	
	private JTextField temperatureField;
	private JTextField NoOfStepsField;
	private JTextField SigmaHSRefFieldA;
	private JTextField SigmaHSRefFieldB;
	
	
	SimulationParameters(SimulationEnvironmentObject ParamObject){
		
			CloseWindow = new JButton("Cancel");
			SaveValues = new JButton("Save");
			GridBagLayout gridbagOtherParam = new GridBagLayout();
			this.setLayout(gridbagOtherParam);
		    
		    temperatureField = new JTextField();
			temperatureField.setColumns(10);
			temperatureField.setSize(new Dimension((int)temperatureField.getSize().getWidth(),20));
			temperatureField.setText(Double.toString(ParamObject.getTemperature()));
			JLabel temperatureLabel = new JLabel("Temperature");
			temperatureLabel.setLabelFor(temperatureField);
			
			Border compoundField1;
			compoundField1 = BorderFactory.createCompoundBorder(
		    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
			temperatureField.setBorder(compoundField1);

			
			
			
			NoOfStepsField = new JTextField();
			NoOfStepsField.setColumns(10);
			JLabel NoOfStepsLabel = new JLabel("Steps");
			NoOfStepsField.setText(Integer.toString(ParamObject.getNoOfSteps()));
			NoOfStepsLabel.setLabelFor(NoOfStepsField);
			
			
			Border compoundField2;
			compoundField2 = BorderFactory.createCompoundBorder(
		    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
			NoOfStepsField.setBorder(compoundField2);
			
			
			SigmaHSRefFieldA = new JTextField();
			SigmaHSRefFieldA.setColumns(10);
			JLabel SigmaHSRefLabelA = new JLabel("SigmaHSRefA");
			SigmaHSRefLabelA.setLabelFor(SigmaHSRefFieldA);
			if(ParamObject.getSigmaHSRefA() == 0.0){
				SigmaHSRefFieldA.setText("0.0");
				SigmaHSRefFieldA.setVisible(false);
				SigmaHSRefLabelA.setVisible(false);
			}else{
			SigmaHSRefFieldA.setText(Double.toString(ParamObject.getSigmaHSRefA()));}
			
			
			Border compoundField3;
			compoundField3 = BorderFactory.createCompoundBorder(
		    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
			SigmaHSRefFieldA.setBorder(compoundField3);
		
			
			SigmaHSRefFieldB = new JTextField();
			SigmaHSRefFieldB.setColumns(10);
			JLabel SigmaHSRefLabelB = new JLabel("SigmaHSRefB");
			SigmaHSRefLabelB.setLabelFor(SigmaHSRefFieldB);
			if(ParamObject.getSigmaHSRefB() == 0.0){
				SigmaHSRefFieldB.setText("0.0");
				SigmaHSRefFieldB.setVisible(false);
				SigmaHSRefLabelB.setVisible(false);
			}else{
			SigmaHSRefFieldB.setText(Double.toString(ParamObject.getSigmaHSRefB()));}
			
			
			Border compoundField4;
			compoundField4 = BorderFactory.createCompoundBorder(
		    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
			SigmaHSRefFieldB.setBorder(compoundField4);

			JComponent[] OtherParamLeft= {temperatureLabel, NoOfStepsLabel, SigmaHSRefLabelA,SigmaHSRefLabelB,SaveValues};
			JComponent[] OtherParamRight = { temperatureField,NoOfStepsField,SigmaHSRefFieldA,SigmaHSRefFieldB,CloseWindow};
		    addLeftRightComponents(OtherParamLeft,OtherParamRight,gridbagOtherParam,this);
		    
		    
		 //   this.add(Message);
		   // this.add(CloseWindow);
		    
		    Border compound2;
		    compound2 = BorderFactory.createCompoundBorder(
		    		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
		    this.setBorder(compound2);

	}
	
	public JTextField getTemperatureField() {
		return temperatureField;
	}

	public JTextField getNoOfStepsField() {
		return NoOfStepsField;
	}

	public JTextField getSigmaHSRefFieldA() {
		return SigmaHSRefFieldA;
	}
	
	public JTextField getSigmaHSRefFieldB() {
		return SigmaHSRefFieldB;
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
		CloseWindow.addActionListener(Comp);
	}

	

	public JButton getCloseWindow() {
		return CloseWindow;
	}

	public JButton getSaveValues() {
		return SaveValues;
	}
	
	

}
