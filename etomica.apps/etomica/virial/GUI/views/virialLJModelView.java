package etomica.virial.GUI.views;

import java.awt.BorderLayout;
import java.awt.ComponentOrientation;
import java.awt.Container;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import java.awt.event.ActionListener;
import javax.swing.*;

import etomica.virial.GUI.models.virialLJModel;

public class virialLJModelView extends JPanel
    {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final String INITIAL_VALUE_nPoints= "3";
	private static final String INITIAL_VALUE_temperature= "2.0";
	private static final String INITIAL_VALUE_steps= "100000";
	private static final String INITIAL_VALUE_sigmaHSRef = "1.5";	
	
	private virialLJModel model;
	//Values for the fields


	//Labels to identify the fields
	private JLabel nPointsLabel;
	private JLabel temperatureLabel;
	private JLabel stepsLabel;
	private JLabel sigmaHSRefLabel;

	//Strings for the labels
	private static String nPointsString = "nPoints: ";
	private static String temperatureString = "Temperature: ";
	private static String stepsString = "Steps: ";
	private static String sigmaHSRefString = "SigmaHSRef: ";

	//Fields for data entry
	private JTextField nPointsField;
	private JTextField temperatureField;
	private JTextField stepsField;
	private JTextField sigmaHSRefField;

	//Buttons for Starting the Simulation and Resetting the values
	private JButton sunSimulation = new JButton("Run");
	private JButton clearValues = new JButton("Reset Values");

	public virialLJModelView(virialLJModel ljModel) {
		super(new BorderLayout());
		model = ljModel;
		model.setValueNPoints(INITIAL_VALUE_nPoints);
		model.setValueTemperature(INITIAL_VALUE_temperature);
		model.setValueSteps(INITIAL_VALUE_steps);
		model.setValueSigmaHSRef(INITIAL_VALUE_sigmaHSRef);

		//Create the labels.
		nPointsLabel = new JLabel(nPointsString);
		
		temperatureLabel = new JLabel(temperatureString);
		 
		stepsLabel = new JLabel(stepsString);
		 
		sigmaHSRefLabel = new JLabel(sigmaHSRefString);
		

		//Create the text fields and set them up.
		nPointsField = new JTextField();
		
		nPointsField.setColumns(10);
		nPointsField.setText(model.getValueNPoints());

		temperatureField = new JTextField();
		temperatureField.setColumns(10);
		temperatureField.setText(model.getValueTemperature());

		stepsField = new JTextField();
		
		stepsField.setColumns(10);
		stepsField.setText(model.getValueSteps());

		sigmaHSRefField = new JTextField();
		 
		sigmaHSRefField.setColumns(10);
		sigmaHSRefField.setText(model.getValueSigmaHSRef());

		nPointsLabel.setLabelFor(nPointsField);
		temperatureLabel.setLabelFor(temperatureField);
		stepsLabel.setLabelFor(stepsField);
		sigmaHSRefLabel.setLabelFor(sigmaHSRefField);

		
		setBorder(BorderFactory.createEmptyBorder(50, 50, 50, 50));
		
		JPanel parametersPanel = new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
        
		parametersPanel.setLayout(gridbag);
        
        JLabel[] labels = {nPointsLabel, temperatureLabel, stepsLabel, sigmaHSRefLabel};
        JTextField[] textFields = {nPointsField, temperatureField, stepsField,sigmaHSRefField};
        addLabelTextRows(labels, textFields, gridbag, parametersPanel);

        JPanel buttonPane = new JPanel();
        FlowLayout buttonPanel = new FlowLayout();
        buttonPane.setLayout(buttonPanel);
        buttonPanel.setAlignment(FlowLayout.TRAILING);
        buttonPane.add(sunSimulation);
        buttonPane.add(clearValues);
        buttonPane.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);

        add(parametersPanel, 
                BorderLayout.PAGE_START);
        add(buttonPane, BorderLayout.LINE_START);
        
	}

	private void addLabelTextRows(JLabel[] labels,
            JTextField[] textFields,
            GridBagLayout gridbag,
            Container container) {
		
			GridBagConstraints c = new GridBagConstraints();
			c.anchor = GridBagConstraints.EAST;
			int numLabels = labels.length;

			for (int i = 0; i < numLabels; i++) {
				c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
				c.fill = GridBagConstraints.NONE;      //reset to default
				c.weightx = 0.0;                       //reset to default
				container.add(labels[i], c);

				c.gridwidth = GridBagConstraints.REMAINDER;     //end row
				c.fill = GridBagConstraints.HORIZONTAL;
				c.weightx = 1.0;
				container.add(textFields[i], c);
			}
}
	//Method to reset the view
	public void reset() {
		nPointsField.setText(INITIAL_VALUE_nPoints);
		temperatureField.setText(INITIAL_VALUE_temperature);
		stepsField.setText(INITIAL_VALUE_steps);
		sigmaHSRefField.setText(INITIAL_VALUE_sigmaHSRef);
	}

	//Method to capture and display the error
	public void showError(String errMessage) {
		JOptionPane.showMessageDialog(this, errMessage);
	}

	//Methods to obtain the user input values for the parameters
	
	public int getUserInputNPoints(){
		return Integer.parseInt(nPointsField.getText());
	}

	public double getUserInputTemperature(){
		return Double.parseDouble(temperatureField.getText());
	}

	public long getUserInputSteps(){
		return Long.parseLong(stepsField.getText());
	}

	public double getUserInputSigmaHSRef(){
		return Double.parseDouble(sigmaHSRefField.getText());
	}
	
	//Listener method to run the simulation when the Run button is pressed
	public void addSendParametersListener(ActionListener mal) {
		sunSimulation.addActionListener(mal);
	}

	//Listener to reset the values when reset Button is pressed
	public void addResetParametersListener(ActionListener cal) {
		clearValues.addActionListener(cal);
	}

}
