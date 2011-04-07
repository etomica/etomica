package etomica.virial.GUI.containers;

import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import etomica.virial.GUI.models.SuperModel;


public class BinaryMixtureParam extends JPanel{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static String componentsA[] = {"LJ","CO2","Methanol","Ethanol","Methane","Ethane","Propane"};
	
	private static String componentsBLJ[] = {"LJ"};
	private static String componentsBCO2Trappe[] = {"Methanol","Ethanol","Methane","Ethane","Propane"};
	private static String componentsBCO2EPM2[] = {"Methane","Ethane","Propane"};
	
	private static String LJPotentialsA[] = {"Sperical-2-Body", "2-Centered-with-Quad"};
	private static String LJPotentialsB[] = {"Sperical-2-Body"};
	
	private static String CO2PotentialsA[] = {"TRAPPE-with-LJ", "2-Centered-with-Quad", "EPM2-with-LJ"};
	private static String CO2PotentialsBTrappe[] = {"TRAPPE-with-LJ"};
	private static String CO2PotentialsBEPM2[] = {"SKS"};
	
	private static String SigmaHSRefLJ = "SigmaHSRef";
	
	
	private static String ComponentA = "Species A";
	private static String ForceFieldA ="Potential/Forcefield";
	
	private static String ComponentB = "Species B";
	private static String ForceFieldB ="Potential/Forcefield";
	
	private static JLabel SigmaHSRefLJLabel;
	private static JLabel ComponentLabelA;
	private static JLabel PotentialLabelA;
	
	private static JLabel ComponentLabelB;
	private static JLabel PotentialLabelB;
	
	private JComboBox componentsListA; 
	private JComboBox componentsListBLJ;
	private JComboBox componentsListBCO2TRAPPE;
	private JComboBox componentsListBCO2EPM2;
	
	
	private JComboBox LJpotentialListA;
	private JComboBox LJpotentialListB;
	private JComboBox CO2PotentialListA;
	
	private JButton DefaultValues;
	private JLabel ButtonText;
	
	private JButton Reset;
	private JLabel ResetLabel;
	
	private JTextField SigmaHSRefLJField;
	private JPanel ComponentBCards;

	public BinaryMixtureParam(){
		super();
		
		GridBagLayout gridbag1 = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		this.setLayout(gridbag1);
		
		ComponentLabelA = new JLabel(ComponentA);
		ComponentLabelB = new JLabel(ComponentB);
		
		componentsListA = new JComboBox(componentsA);
		
		ComponentBCards = new JPanel(new CardLayout());
		componentsListBLJ = new JComboBox(componentsBLJ);
		componentsListBCO2TRAPPE = new JComboBox(componentsBCO2Trappe);
		componentsListBCO2EPM2 = new JComboBox(componentsBCO2EPM2);
		
		ComponentBCards.add("CompBLJ",componentsListBLJ);
		ComponentBCards.add("CompBCO2Trappe",componentsListBCO2TRAPPE);
		ComponentBCards.add("CompBCO2EMP2",componentsListBCO2EPM2);
		
		
		
		
		JComponent[] components1 = {ComponentLabelA,componentsListA};
		JComponent[] components2 = {ComponentLabelB,ComponentBCards};
		addLabelWithComponents(components1,components2,gridbag1,this);
	}
	
	private void addLabelWithComponents(JComponent[] components1,
            JComponent[] components2,
            GridBagLayout gridbag,
            Container container) {

			GridBagConstraints c = new GridBagConstraints();
			//c.anchor = GridBagConstraints.WEST;
			
			
			int numLabels = components1.length;

			for (int i = 0; i < numLabels; i++) {
		
			c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
			//c.fill = GridBagConstraints.NONE;      //reset to default
			c.fill = GridBagConstraints.BOTH;
			c.insets = new Insets(10,10,10,10);
			c.ipady = 5;
			c.anchor = GridBagConstraints.NORTHEAST;
			c.weightx = 0.5;     
			
			
			//reset to default
			container.add(components1[i], c);

			c.gridwidth = GridBagConstraints.REMAINDER;     //end row
			c.fill = GridBagConstraints.BOTH;
			c.anchor = GridBagConstraints.EAST;
			//c.weightx = 1.0;
			container.add(components2[i], c);}
		}

	public void addChooseComponentAListener(ActionListener Comp) {
		componentsListA.addActionListener(Comp);
	}

	public JPanel getComponentBCards() {
		return ComponentBCards;
	}

	public void setComponentBCards(JPanel componentBCards) {
		ComponentBCards = componentBCards;
	}
}
