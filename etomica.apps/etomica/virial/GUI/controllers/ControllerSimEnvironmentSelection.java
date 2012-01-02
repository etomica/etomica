package etomica.virial.GUI.controllers;


import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import etomica.virial.GUI.models.ModelSimulationEnvironment;
import etomica.virial.GUI.models.ModelTemperatureAndSteps;
import etomica.virial.GUI.views.ViewSimEnvironmentSelection;


public class ControllerSimEnvironmentSelection {
	private ViewSimEnvironmentSelection simEnvironmentView;
	
	private ModelSimulationEnvironment simEnvironmentDM;
	private ModelTemperatureAndSteps modelTemperatureAndSteps;
	
	
	public ControllerSimEnvironmentSelection(ModelSimulationEnvironment simEnvObject){
		this.simEnvironmentDM = simEnvObject;
		modelTemperatureAndSteps = ModelTemperatureAndSteps.getInstance();
		simEnvironmentView = new ViewSimEnvironmentSelection(this.simEnvironmentDM);
		simEnvironmentView.openSimEnvParamView();
		simEnvironmentView.addSaveValuesButtonListener(new SaveValuesButtonListener());
		simEnvironmentView.addCloseWindowButtonListener(new CloseWindowButtonListener());
	}
	
	
	class SaveValuesButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			modelTemperatureAndSteps.setTemperature(Double.parseDouble(simEnvironmentView.getTemperatureField().getText()));
			simEnvironmentDM.setTemperature(Double.parseDouble(simEnvironmentView.getTemperatureField().getText()));
			modelTemperatureAndSteps.setNoOfSteps(Integer.parseInt(simEnvironmentView.getNoOfStepsField().getText()));
			simEnvironmentDM.setNoOfSteps(modelTemperatureAndSteps.getNoOfSteps());
			for(int i=0;i<10;i++){
				simEnvironmentDM.setSigmaHSRef(Double.parseDouble(simEnvironmentView.getSigmaHSRefField(i).getText()),i);
			}
			simEnvironmentView.closeSimEnvParamView();
		
		}

	
	}
	
	class CloseWindowButtonListener implements ActionListener 
	 {

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			modelTemperatureAndSteps.reset();
			simEnvironmentView.closeSimEnvParamView();
		}

	
	}
	
}
