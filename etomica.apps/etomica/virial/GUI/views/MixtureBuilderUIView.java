package etomica.virial.GUI.views;

import etomica.virial.GUI.containers.AlertMessageUIView;
import etomica.virial.GUI.containers.AlkaneSpheresUIView;
import etomica.virial.GUI.containers.MixtureBuilderUIFrame;
import etomica.virial.GUI.containers.SimulationEnvironmentUIView;
import etomica.virial.GUI.containers.SpeciesSelectionUIView;
import etomica.virial.GUI.models.MixtureBuilderDM;


public class MixtureBuilderUIView {
		
	
	private SpeciesSelectionUIView speciesSelectionView;
	private AlkaneSpheresUIView alkaneLengthSelectionView1;
	private AlkaneSpheresUIView alkaneLengthSelectionView2;
	private AlertMessageUIView alertMessageView;
	private SimulationEnvironmentUIView simEnvView;
	
	
	
	public MixtureBuilderUIView(MixtureBuilderUIFrame mixtureBuilder, MixtureBuilderDM SM) {
		
		
		
		speciesSelectionView = new SpeciesSelectionUIView();
	    mixtureBuilder.add(speciesSelectionView);
	    mixtureBuilder.setFrameProperties();
	    reset();
	}


	public void reset(){
		
	}



}
