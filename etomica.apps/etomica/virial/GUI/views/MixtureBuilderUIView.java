package etomica.virial.GUI.views;

import etomica.virial.GUI.containers.AlertMessageUIView;
import etomica.virial.GUI.containers.MixtureBuilderUIFrame;
import etomica.virial.GUI.containers.SpeciesSelectionUIView;
import etomica.virial.GUI.models.MixtureBuilderDM;


public class MixtureBuilderUIView {
		
	
	private SpeciesSelectionUIView speciesSelectionView;
	private AlertMessageUIView DialogBox;
	
	
	
	public MixtureBuilderUIView(MixtureBuilderUIFrame mixtureBuilder, MixtureBuilderDM SM) {
		
		
		
		speciesSelectionView = new SpeciesSelectionUIView();
	    mixtureBuilder.add(speciesSelectionView);
	    mixtureBuilder.setFrameProperties();
	    reset();
	}


	public void reset(){
		
	}


	public AlertMessageUIView getDialogBox() {
		return DialogBox;
	}

	public void setDialogBox(AlertMessageUIView dialogBox) {
		DialogBox = dialogBox;
	}

	
	public SpeciesSelectionUIView getSpeciesParameters() {
		return speciesSelectionView;
	}

	public void setSpeciesParameters(SpeciesSelectionUIView SpeciesSelectionView) {
		speciesSelectionView = SpeciesSelectionView;
	}
	


}
