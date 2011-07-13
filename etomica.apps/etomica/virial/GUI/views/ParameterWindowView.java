package etomica.virial.GUI.views;

import etomica.virial.GUI.containers.DialogBoxPanel;
import etomica.virial.GUI.containers.MainFrame;
import etomica.virial.GUI.containers.Species1;
import etomica.virial.GUI.containers.SuperContainer;
import etomica.virial.GUI.models.SuperModel;


public class ParameterWindowView {
		
	
	private Species1 species1;
	private DialogBoxPanel DialogBox;
	
	
	
	public ParameterWindowView(MainFrame mixtureBuilder, SuperModel SM) {
		
		SuperContainer container = new SuperContainer();
		
		species1 = container.getS1();
	    mixtureBuilder.add(species1);

	    mixtureBuilder.SetFrameProperties();
	    reset();
	}


	public void reset(){
		
	}


	public DialogBoxPanel getDialogBox() {
		return DialogBox;
	}

	public void setDialogBox(DialogBoxPanel dialogBox) {
		DialogBox = dialogBox;
	}

	
	public Species1 getSpeciesParameters() {
		return species1;
	}

	public void setSpeciesParameters(Species1 speciesParameters) {
		species1 = speciesParameters;
	}
	


}
