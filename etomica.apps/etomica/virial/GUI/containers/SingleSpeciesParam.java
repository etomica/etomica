package etomica.virial.GUI.containers;

import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionListener;
import java.awt.event.ItemListener;
import java.text.DecimalFormat;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.SoftBevelBorder;

import etomica.virial.GUI.models.SingleSpeciesModel;
import etomica.virial.GUI.models.SuperModel;


public class SingleSpeciesParam extends JPanel{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	

	private static String components[] = {"LJ","CO2","Methanol","Ethanol","Methane","Ethane","Propane"};
	private static String LJPotentials[] = {"Spherical-2-Body", "2-Centered-with-Quad", "Spherical-2-Body-with-Quad"};
	private static String CO2Potentials[] = {"TRAPPE-with-LJ", "2-Centered-with-Quad", "EPM2-with-LJ"};
	private static String AlkanePotentials[] = {"TRAPPE-with-LJ", "SKS-with-LJ", "EPM2-with-LJ"};
	private static String AlcoholPotentials[] = {"TRAPPE-with-LJ", "ROWLEY"};
	
	
	private static String SigmaHSRefLJ = "SigmaHSRef";
	
	
	private static String Component = "Species A";
	private static String ForceField ="Potential/Forcefield";
	
	private static JLabel ComponentLabel;
	private static JLabel PotentialLabel;
	private static JLabel SigmaHSRefLJLabel;
	
	public static JLabel getPotentialLabel() {
		return PotentialLabel;
	}

	public static void setPotentialLabel(JLabel potentialLabel) {
		PotentialLabel = potentialLabel;
	}

	private JComboBox componentsList; 
	private JComboBox LJpotentialList;
	private JComboBox CO2PotentialList;
	private JComboBox AlkanePotentialList;
	private JComboBox AlcoholPotentialsList;
	
	private JButton DefaultValues;
	private JLabel ButtonText;
	
	private JButton Reset;
	private JLabel ResetLabel;
	
	private JTextField SigmaHSRefLJField;
	
	
	//Card Layout Panels
	private JPanel PotentialCards;
	private JPanel SigmaParametersLJ;
	
	
	
	public SingleSpeciesParam(){
		super();
		
		GridBagLayout gridbag1 = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		this.setLayout(gridbag1);
		
		ComponentLabel = new JLabel(Component);
		componentsList = new JComboBox(components);
		componentsList.setToolTipText("Select one of the components on the dropdown");
		//componentsList.setSelectedIndex(0);
		
		
		PotentialLabel = new JLabel(ForceField);
		
		
		
		LJpotentialList = new JComboBox(LJPotentials);
		CO2PotentialList = new JComboBox(CO2Potentials);
		AlkanePotentialList = new JComboBox(AlkanePotentials);
		AlcoholPotentialsList = new JComboBox(AlcoholPotentials);
		
		
		PotentialCards = new JPanel(new CardLayout());
		
		PotentialCards.add("LJ",LJpotentialList);
		PotentialCards.add("CO2",CO2PotentialList);
		PotentialCards.add("Alkane",AlkanePotentialList);
		PotentialCards.add("Alcohol",AlcoholPotentialsList);
		
		
		SigmaHSRefLJLabel = new JLabel(SigmaHSRefLJ);
		SigmaHSRefLJField = new JTextField(10);
		
		ButtonText = new JLabel("Press to edit");
		DefaultValues = new JButton("Default Values");
		
		ResetLabel = new JLabel("Press to reset Values");
		Reset = new JButton("Reset");
		
		this.add(DefaultValues);
		
		
		
		JLabel[] labels = {ComponentLabel,PotentialLabel,SigmaHSRefLJLabel,ButtonText,ResetLabel};
		JComponent[] components = {componentsList,PotentialCards,SigmaHSRefLJField,DefaultValues,Reset};
		addLabelWithComponents(labels,components,gridbag1, this);
		
		c.gridwidth = GridBagConstraints.REMAINDER; //last
	    c.anchor = GridBagConstraints.WEST;
	}
	
	public JButton getReset() {
		return Reset;
	}

	public void setReset(JButton reset) {
		Reset = reset;
	}

	public JLabel getButtonText() {
		return ButtonText;
	}

	public JButton getDefaultValues() {
		return DefaultValues;
	}

	public void setDefaultValues(JButton defaultValues) {
		DefaultValues = defaultValues;
	}

	public static JLabel getSigmaHSRefLJLabel() {
		return SigmaHSRefLJLabel;
	}

	public static void setSigmaHSRefLJLabel(JLabel sigmaHSRefLJLabel) {
		SigmaHSRefLJLabel = sigmaHSRefLJLabel;
	}

	public JTextField getSigmaHSRefLJField() {
		return SigmaHSRefLJField;
	}

	public void setSigmaHSRefLJField(JTextField sigmaHSRefLJField) {
		SigmaHSRefLJField = sigmaHSRefLJField;
	}

	private void addLabelWithComponents(JLabel[] labels,
            JComponent[] components,
            GridBagLayout gridbag,
            Container container) {

			GridBagConstraints c = new GridBagConstraints();
			//c.anchor = GridBagConstraints.WEST;
			
			
			int numLabels = labels.length;

			for (int i = 0; i < numLabels; i++) {
		
			c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
			//c.fill = GridBagConstraints.NONE;      //reset to default
			c.fill = GridBagConstraints.BOTH;
			c.insets = new Insets(10,10,10,10);
			c.ipady = 5;
			c.anchor = GridBagConstraints.NORTHEAST;
			c.weightx = 0.5;     
			
			
			//reset to default
			container.add(labels[i], c);

			c.gridwidth = GridBagConstraints.REMAINDER;     //end row
			c.fill = GridBagConstraints.BOTH;
			c.anchor = GridBagConstraints.EAST;
			//c.weightx = 1.0;
			container.add(components[i], c);}
		}

	
	
	
	public JPanel getPotentialCards() {
		return PotentialCards;
	}

	public void setPotentialCards(JPanel potentialCards) {
		PotentialCards = potentialCards;
	}

	int getUserChoiceComponent(){
		return componentsList.getSelectedIndex();
	}
	
	int getUserChoiceLJPotential(){
		return LJpotentialList.getSelectedIndex();
	}
	
	int getUserChoiceCO2Potential(){
		return CO2PotentialList.getSelectedIndex();
	}
	
	public void setSigmaHSRefLJText(double SigmaHSRefLJ){
		SigmaHSRefLJField.setText(Double.toString(SigmaHSRefLJ));
	}
	
	
	public void addChooseComponentListener(ActionListener Comp) {
		componentsList.addActionListener(Comp);
	}
	
	public void addChooseLJListener(ActionListener Comp) {
		LJpotentialList.addActionListener(Comp);
	}
	
	public void addChooseCO2Listener(ActionListener Comp) {
		CO2PotentialList.addActionListener(Comp);
	}
	
	public void addChooseAlkaneListener(ActionListener Comp) {
		AlkanePotentialList.addActionListener(Comp);
	}
	
	public void addChooseAlcoholListener(ActionListener Comp) {
		AlcoholPotentialsList.addActionListener(Comp);
	}
	
	public void addChooseResetListener(ActionListener Comp) {
		Reset.addActionListener(Comp);
	}
	
	public void addChooseDefaultListener(ActionListener Comp) {
		DefaultValues.addActionListener(Comp);
	}
	
	public JLabel getResetLabel() {
		return ResetLabel;
	}

	public void setResetLabel(JLabel resetLabel) {
		ResetLabel = resetLabel;
	}

	public JComboBox getComponentsList() {
		return componentsList;
	}

	public void setComponentsList(JComboBox componentsList) {
		this.componentsList = componentsList;
	}

	public JComboBox getLJpotentialList() {
		return LJpotentialList;
	}

	public void setLJpotentialList(JComboBox lJpotentialList) {
		LJpotentialList = lJpotentialList;
	}

	public JComboBox getCO2PotentialList() {
		return CO2PotentialList;
	}

	public void setCO2PotentialList(JComboBox cO2PotentialList) {
		CO2PotentialList = cO2PotentialList;
	}

	public JComboBox getAlkanePotentialList() {
		return AlkanePotentialList;
	}

	public void setAlkanePotentialList(JComboBox alkanePotentialList) {
		AlkanePotentialList = alkanePotentialList;
	}

	public JComboBox getAlcoholPotentialsList() {
		return AlcoholPotentialsList;
	}

	public void setAlcoholPotentialsList(JComboBox alcoholPotentialsList) {
		AlcoholPotentialsList = alcoholPotentialsList;
	}
	
}
