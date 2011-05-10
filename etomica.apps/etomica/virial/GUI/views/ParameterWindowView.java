package etomica.virial.GUI.views;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionListener;

import javax.swing.JFrame;


import etomica.virial.GUI.containers.BinaryMixtureParam;
import etomica.virial.GUI.containers.DialogBoxPanel;
import etomica.virial.GUI.containers.LJDefaultViewPanel;
import etomica.virial.GUI.containers.MainFrame;
import etomica.virial.GUI.containers.MainFramePanel;
import etomica.virial.GUI.containers.ParameterTabs;
import etomica.virial.GUI.containers.RunParam;
import etomica.virial.GUI.containers.SingleSpeciesParam;
import etomica.virial.GUI.containers.Species;
import etomica.virial.GUI.containers.Species1;
import etomica.virial.GUI.containers.SuperContainer;
import etomica.virial.GUI.models.RunParametersModel;
import etomica.virial.GUI.models.SingleSpeciesModel;
import etomica.virial.GUI.models.SpeciesModel;
import etomica.virial.GUI.models.SuperModel;
import etomica.virial.simulations.VirialAlkane;
import etomica.virial.simulations.VirialLJ;


public class ParameterWindowView {
		
	private ParameterTabs ParamTabs;
	private Species1 SpeciesParameters;
	private RunParam RunParameters;
	private DialogBoxPanel DialogBox;
	private LJDefaultViewPanel LJDefault;
	private JFrame DialogFrame;
	private JFrame DefaultViewFrame;
	
	
	public ParameterWindowView(MainFrame mixtureBuilder, SuperModel SM) {
		
		SuperContainer container = new SuperContainer();
		
		DialogBox = container.getDBox();
		DialogFrame = new JFrame("Message");
		DialogFrame.add(DialogBox);
		DialogFrame.setVisible(false);
		
		DefaultViewFrame = new JFrame("Default Values");
		LJDefault = container.getLJD();
		DefaultViewFrame.add(LJDefault);
		DefaultViewFrame.setVisible(false);
		
		
		
		//Instantiate the Main Frame and properties.
		mixtureBuilder.setLayout(new GridLayout(1,3));
		
		//Instantiate the MainFrame Panels - Parameters , Graphic and Console
		MainFramePanel ParameterPanel = new MainFramePanel(900, 800);
		mixtureBuilder.add(ParameterPanel);
		
		SpeciesParameters = container.getS1();
		RunParameters = container.getRP();
	
		ParamTabs = new ParameterTabs(4);
		ParamTabs.MakeIndividualTabs(ParamTabs.getTabs(),4,SpeciesParameters,RunParameters);
		
	    ParameterPanel.add(ParamTabs);
	    
	    
	    MainFramePanel GraphicPanel = new MainFramePanel(280, 800);
	    mixtureBuilder.add(GraphicPanel);

	    mixtureBuilder.SetFrameProperties();
	    reset();
	}

	public LJDefaultViewPanel getLJDefault() {
		return LJDefault;
	}

	public void setLJDefault(LJDefaultViewPanel lJDefault) {
		LJDefault = lJDefault;
	}

	public JFrame getDefaultViewFrame() {
		return DefaultViewFrame;
	}

	public void setDefaultViewFrame(JFrame defaultViewFrame) {
		DefaultViewFrame = defaultViewFrame;
	}

	public void reset(){
		
	}

	public JFrame getDialogFrame() {
		return DialogFrame;
	}

	public void setDialogFrame(JFrame dialogFrame) {
		DialogFrame = dialogFrame;
	}

//DialogBoxPanel GEtters and Setters


	public DialogBoxPanel getDialogBox() {
		return DialogBox;
	}

	public void setDialogBox(DialogBoxPanel dialogBox) {
		DialogBox = dialogBox;
	}

//ParamTabs Getter and Setters

	public ParameterTabs getParamTabs() {
		return ParamTabs;
	}

	public void setParamTabs(ParameterTabs paramTabs) {
		ParamTabs = paramTabs;
	}

	
	public Species1 getSpeciesParameters() {
		return SpeciesParameters;
	}

	public void setSpeciesParameters(Species1 speciesParameters) {
		SpeciesParameters = speciesParameters;
	}
	

	public RunParam getRunParameters() {
		return RunParameters;
	}

	public void setRunParameters(RunParam runParameters) {
		RunParameters = runParameters;
	}


	public ParameterTabs getParameterPanelTabs() {
		return ParamTabs;
	}

	public void setParameterPanelTabs(ParameterTabs parameterPanelTabs) {
		ParamTabs = parameterPanelTabs;
	}


}
