package etomica.virial.GUI;

import javax.swing.SwingUtilities;

import etomica.virial.GUI.containers.MainFrame;
import etomica.virial.GUI.controllers.ParameterController;
import etomica.virial.GUI.models.RunParametersModel;
import etomica.virial.GUI.models.SingleSpeciesModel;
import etomica.virial.GUI.models.SpeciesModel;
import etomica.virial.GUI.models.SuperModel;
import etomica.virial.GUI.views.ParameterWindowView;
import etomica.virial.GUI.containers.MainFramePanel;



public class MixtureBuilder{
	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				MainFrame MixtureBuilder = new MainFrame("Mixture Builder");
				//MainFrame DefaultValueFrame = new MainFrame("Default Values");
			
				SuperModel superModel = new SuperModel();
				//Instantiate the view for MixtureParameters
				//DefaultValuesView DWindow = new DefaultValuesView(DefaultValueFrame,superModel);
				ParameterWindowView PWindow = new ParameterWindowView(MixtureBuilder,superModel);
				
				@SuppressWarnings("unused")
				ParameterController controller = new ParameterController(PWindow,superModel);
				
				System.out.println("Hi THIS IS CONNECTED");
				
			}
		});
    }

}
