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
import etomica.virial.GUI.containers.SuperContainer;
import etomica.virial.GUI.models.RunParametersModel;
import etomica.virial.GUI.models.SingleSpeciesModel;
import etomica.virial.GUI.models.SpeciesModel;
import etomica.virial.GUI.models.SuperModel;
import etomica.virial.simulations.VirialAlkane;
import etomica.virial.simulations.VirialLJ;


public class ParameterWindowView {
		
	private ParameterTabs ParamTabs;
	private Species SpeciesParameters;
	private SingleSpeciesParam OneSpecies;
	private BinaryMixtureParam TwoSpecies;
	private RunParam RunParameters;
	private DialogBoxPanel DialogBox;
	private LJDefaultViewPanel LJDefault;
	private JFrame DialogFrame;
	private JFrame DefaultViewFrame;
	private SpeciesModel speciesmodel;
	private SingleSpeciesModel onespeciesmodel;
	
	public ParameterWindowView(MainFrame mixtureBuilder, SuperModel SM) {
		
		SuperContainer container = new SuperContainer();
		speciesmodel = SM.getSm();
		onespeciesmodel = SM.getSsm();
		speciesmodel.setPureBinaryButtonValue(SpeciesModel.getIntialValuePureorbinaryrbuttonvalue());
		
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
		MainFramePanel ParameterPanel = new MainFramePanel(510, 616);
		mixtureBuilder.add(ParameterPanel);
		
		
		OneSpecies = container.getSSP();
		
		
		TwoSpecies = container.getBMP();
		
		
		SpeciesParameters = container.getS();
		SpeciesParameters.setOneSpecies(OneSpecies);
		SpeciesParameters.setTwoSpecies(TwoSpecies);
		SpeciesParameters.addCards(OneSpecies, TwoSpecies);
		
		RunParameters = container.getRP();
	
		ParamTabs = new ParameterTabs(4, 510);
		ParamTabs.MakeIndividualTabs(ParamTabs.getTabs(),4,SpeciesParameters,RunParameters);
		
		if(speciesmodel.getPureBinaryButtonValue()==0){
			SpeciesParameters.getPure().setSelected(true);}
		if(speciesmodel.getPureBinaryButtonValue()==1){
			SpeciesParameters.getBinary().setSelected(true);}
	    ParameterPanel.add(ParamTabs);
	    
	    
	    MainFramePanel GraphicPanel = new MainFramePanel(280, 616);
	    mixtureBuilder.add(GraphicPanel);

	    mixtureBuilder.SetFrameProperties();
	  //  VirialAlkane.VirialSiepmannSpheresParam lj = new VirialAlkane.VirialSiepmannSpheresParam(3,3,300.0,1000,-1,false);
	   // VirialAlkane.runVirialAlkane(lj, new JFrame());
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
		OneSpecies.getPotentialCards().setVisible(false);
		SingleSpeciesParam.getPotentialLabel().setVisible(false);
		OneSpecies.setSigmaHSRefLJText(SingleSpeciesModel.getInitialValueSigmahsreflj());
		OneSpecies.getSigmaHSRefLJField().setVisible(false);
		SingleSpeciesParam.getSigmaHSRefLJLabel().setVisible(false);
		OneSpecies.getDefaultValues().setVisible(false);
		OneSpecies.getButtonText().setVisible(false);
		OneSpecies.getResetLabel().setVisible(false);
		OneSpecies.getReset().setVisible(false);
		TwoSpecies.getComponentBCards().setVisible(false);
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


	public Species getSpeciesParameters() {
		return SpeciesParameters;
	}

	public void setSpeciesParameters(Species speciesParameters) {
		SpeciesParameters = speciesParameters;
	}


	public SingleSpeciesParam getOneSpecies() {
		return OneSpecies;
	}

	public void setOneSpecies(SingleSpeciesParam oneSpecies) {
		OneSpecies = oneSpecies;
	}

	public BinaryMixtureParam getTwoSpecies() {
		return TwoSpecies;
	}

	public void setTwoSpecies(BinaryMixtureParam twoSpecies) {
		TwoSpecies = twoSpecies;
	}

	public RunParam getRunParameters() {
		return RunParameters;
	}

	public void setRunParameters(RunParam runParameters) {
		RunParameters = runParameters;
	}


	public SpeciesModel getSpeciesmodel() {
		return speciesmodel;
	}

	public void setSpeciesmodel(SpeciesModel speciesmodel) {
		this.speciesmodel = speciesmodel;
	}

	public ParameterTabs getParameterPanelTabs() {
		return ParamTabs;
	}

	public void setParameterPanelTabs(ParameterTabs parameterPanelTabs) {
		ParamTabs = parameterPanelTabs;
	}


}
