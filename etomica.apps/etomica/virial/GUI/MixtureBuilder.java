package etomica.virial.GUI;

import javax.swing.SwingUtilities;

import etomica.virial.GUI.containers.MixtureBuilderUIFrame;

import etomica.virial.GUI.models.MixtureBuilderDM;
import etomica.virial.GUI.views.MixtureBuilderUIView;




public class MixtureBuilder{
	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				MixtureBuilderUIFrame mixtureBuilderFrame = new MixtureBuilderUIFrame("Mixture Builder");
				//MainFrame DefaultValueFrame = new MainFrame("Default Values");
			
				MixtureBuilderDM mixtureBuilderDM = new MixtureBuilderDM();
				//Instantiate the view for MixtureParameters
				//DefaultValuesView DWindow = new DefaultValuesView(DefaultValueFrame,superModel);
				MixtureBuilderUIView mixtureBuilderView = new MixtureBuilderUIView(mixtureBuilderFrame,mixtureBuilderDM);
				
				
				
				
				System.out.println("Hi THIS IS CONNECTED");
				
			}
		});
    }

}
