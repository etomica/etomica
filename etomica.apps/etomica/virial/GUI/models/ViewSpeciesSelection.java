/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;


import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerListModel;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;

import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;


import javax.swing.border.Border;


import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;
import javax.swing.text.TableView.TableRow;







import etomica.virial.GUI.models.ModelMoleculeListDisplayTable;
import etomica.virial.GUI.models.ModelMoleculeListDisplayTable_new;
import etomica.virial.GUI.models.ModelSelectedSpecies;
import etomica.virial.GUI.models.ModelPotParamDisplayTable;
import etomica.virial.GUI.models.ModelSpeciesSelection;


public class ViewSpeciesSelection extends JPanel {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	//All DataModels this view requires
	private ModelSpeciesSelection modelSpeciesSelection;
	private ModelSelectedSpecies modelSelectedSpecies;
	
	//All JButtons seen on this view
	private JButton AddSpecies;
	private JButton AddPotential;
	private JButton EditVariables;
	
	private JButton AddSameSpecies;
	private JButton AddAnotherSpecies;
	private JButton Remove;
	
	private JButton Reset;
	private JButton AlterSimEnv;
	private JButton Run;
	
	private JButton addSpecies1;
	private JButton addSpecies2;
	private JButton addSpecies3;
	
	private JButton removeSpecies1;
	private JButton removeSpecies2;
	private JButton removeSpecies3;
	
	private JTextField speciesCount1;
	private JTextField speciesCount2;
	private JTextField speciesCount3;

	//Parameter Tables
	private JTable potentialParamDispayTable;
	private JTable moleculeListDisplayTable;
	
	private ListSelectionModel speciesListDisplayLSM;
	private ListSelectionModel potentialListDisplayLSM;
	
	private ModelPotParamDisplayTable potentialParamDisplayTM;
	private TableModelListener potentialParamDisplayTMListener;
	
	private ModelMoleculeListDisplayTable_new moleculeListDisplayTM;
	private TableModelListener moleculeListDisplayTMListener;
	
	private ListSelectionModel potentialParamDisplayLSM;
	private ListSelectionListener potentialParamDisplayLSMListener;
	
	private ListSelectionModel moleculeListDisplayLSM;
	private ListSelectionListener moleculeListDisplayLSMListener;
	
	
	//Description
	private JTextArea Description; 
	//JList 
	private JList speciesJList;
	private ListSelectionListener speciesJListListener;
	
	private JList potentialJList;
	private ListSelectionListener potentialJListListener;
	
	public ViewSpeciesSelection(ModelSpeciesSelection modelSpeciesSelection, ModelSelectedSpecies modelSelectedSpecies){
		super();
		this.modelSpeciesSelection = modelSpeciesSelection;
		this.modelSelectedSpecies = modelSelectedSpecies;
		initComponents();
		}
		
		public void initComponents(){
		
		this.setLayout(new BorderLayout());
	
		JPanel panel1 = new JPanel();
		GridBagLayout gridbaglayout = new GridBagLayout();
		panel1.setLayout(gridbaglayout);
		
		this.add(panel1,BorderLayout.NORTH);
	
		AddSpecies = new JButton("Add");
		AddSpecies.setFocusPainted(false);
		
		speciesJList = new JList(ModelSpeciesSelection.getIntialdisplaylist());
		speciesJList.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		speciesListDisplayLSM = speciesJList.getSelectionModel();
		
		speciesJList.setLayoutOrientation(JList.HORIZONTAL_WRAP);
		speciesJList.setVisibleRowCount(-1);
		speciesJList.setBorder(BorderFactory.createLineBorder(Color.black));
	
 
		JScrollPane SpeciesListScroller = new JScrollPane(speciesJList);
		SpeciesListScroller.setPreferredSize(new Dimension(400,100));
		SpeciesListScroller.setBorder(BorderFactory.createLoweredBevelBorder());
    
		JPanel listPane = new JPanel();				
		listPane.setLayout(new BoxLayout(listPane, BoxLayout.PAGE_AXIS));
    
		JLabel label = new JLabel("Species");
		label.setForeground(Color.WHITE);
		JPanel labelPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		labelPanel.setBackground(new Color(0,0,225));
		labelPanel.setBounds(listPane.getBounds().x, listPane.getBounds().y, 400, label.getHeight());

		label.setLabelFor(speciesJList);
		labelPanel.add(label);
    
		listPane.add(labelPanel);
		listPane.add(Box.createRigidArea(new Dimension(0,0)));
		listPane.add(SpeciesListScroller);
    
    
		Border compound1;
		compound1 = BorderFactory.createCompoundBorder(
				BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
		listPane.setBorder(compound1);
    

		AddPotential = new JButton("Add");
		AddPotential.setVisible(true);
		AddPotential.setEnabled(false);
		AddPotential.setFocusPainted(false);
   
		potentialJList = new JList(ModelSpeciesSelection.getIntialpotentialdisplaylist());
		potentialJList.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		potentialJList.setLayoutOrientation(JList.HORIZONTAL_WRAP);
		potentialJList.setVisibleRowCount(-1);
		potentialJList.setEnabled(false);
 
		JScrollPane PotentialListScroller = new JScrollPane(potentialJList);
		PotentialListScroller.setPreferredSize(new Dimension(400,100));
		PotentialListScroller.setBorder(BorderFactory.createLoweredBevelBorder());
  
		JPanel list2Pane = new JPanel();
		list2Pane.setLayout(new BoxLayout(list2Pane, BoxLayout.PAGE_AXIS));
   
		JLabel label2 = new JLabel("Potentials");
		label2.setForeground(Color.WHITE);
		label2.setLabelFor(potentialJList);
    
		JPanel labelPanel2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		labelPanel2.setBackground(new Color(0,0,225));
		labelPanel2.setBounds(list2Pane.getBounds().x, list2Pane.getBounds().y, 400, label2.getHeight());
		labelPanel2.add(label2);
    
		list2Pane.add(labelPanel2);
		list2Pane.add(Box.createRigidArea(new Dimension(5,0)));
		list2Pane.add(PotentialListScroller);
    
		Border compound2;
		compound2 = BorderFactory.createCompoundBorder(
				BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
		list2Pane.setBorder(compound2);
    
		EditVariables = new JButton("Edit");
		EditVariables.setVisible(true);
		EditVariables.setEnabled(false);
		EditVariables.setFocusPainted(false);
		//EditVariables.addActionListener(this);
    
		JPanel ParametersPane1 = new JPanel();
		ParametersPane1.setLayout(new BoxLayout(ParametersPane1, BoxLayout.PAGE_AXIS));
		
		JPanel ParametersPane2 = new JPanel();
		ParametersPane2.setLayout(new BoxLayout(ParametersPane2, BoxLayout.PAGE_AXIS));
		
		JPanel ParametersPane = new JPanel();
		ParametersPane.setLayout(new GridLayout(0,2));
		ParametersPane.add(ParametersPane1);
		ParametersPane.add(ParametersPane2);
    
    
		JPanel TablePane1 = new JPanel();
		TablePane1.setBorder(BorderFactory.createLoweredBevelBorder());
    
		JPanel TablePane2 = new JPanel();
		TablePane2.setBorder(BorderFactory.createLoweredBevelBorder());
    
		potentialParamDisplayTM = new ModelPotParamDisplayTable();
		moleculeListDisplayTM = new ModelMoleculeListDisplayTable_new();
		

		potentialParamDispayTable = new JTable();
		potentialParamDispayTable.setModel(potentialParamDisplayTM);
		potentialParamDispayTable.setEnabled(false);
		potentialParamDisplayLSM = potentialParamDispayTable.getSelectionModel();
		
		
		
		
		TableColumn column1 = null;
		for (int i = 0; i < 2; i++) {
			column1 = potentialParamDispayTable.getColumnModel().getColumn(i);
			
			column1.setPreferredWidth(150);
		}
    
        moleculeListDisplayTable = new JTable();
        
        moleculeListDisplayTable.setModel(moleculeListDisplayTM);
        moleculeListDisplayTable.setEnabled(false);
		moleculeListDisplayLSM = moleculeListDisplayTable.getSelectionModel();
		
        
    
        
        TableColumn column2 = null;
        for (int i = 0; i < 2; i++) {
        	column2 = moleculeListDisplayTable.getColumnModel().getColumn(i);
        	if(i==0){
        		column2.setPreferredWidth(80);}
        	else{
        		column2.setPreferredWidth(150);
        	}
       
        }
        moleculeListDisplayTable.setRowHeight(28);
        /*moleculeListDisplayTable.setRowMargin(100);
        System.out.println(moleculeListDisplayTable.getIntercellSpacing().width);
        System.out.println(moleculeListDisplayTable.getIntercellSpacing().height);
        moleculeListDisplayTable.setGridColor(Color.red);*/
        TablePane1.add(potentialParamDispayTable);
        TablePane2.add(moleculeListDisplayTable);
        
        
        
    
        JLabel label3 = new JLabel("Parameters");
        label3.setForeground(Color.WHITE);
        JPanel labelPanel3 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        labelPanel3.setMaximumSize(new Dimension(500,label3.getHeight()));
        //labelPanel3.setBackground(new Color(0,78,152));
        labelPanel3.setBackground(new Color(0,0,225));
    
        /*labelPanel3.setBounds(ParametersPane.getBounds().x, ParametersPane.getBounds().y, 400, label3.getHeight());*/
        labelPanel3.setBounds(ParametersPane1.getBounds().x, ParametersPane2.getBounds().y, 400, label3.getHeight());
        labelPanel3.add(label3);
    
        JLabel label4 = new JLabel("Species List");
        label4.setForeground(Color.WHITE);
        
        JPanel labelPanel4 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        labelPanel4.setMaximumSize(new Dimension(500,label4.getHeight()));
        labelPanel4.setBackground(new Color(0,0,225));
        
        labelPanel4.setBounds(ParametersPane.getBounds().x + labelPanel3.getWidth(), ParametersPane.getBounds().y,500, label4.getHeight());
        labelPanel4.add(label4);
    
        //JPanel labelPane = new JPanel();
        //labelPane.add(labelPanel3);
        //labelPane.add(labelPanel4);
   
    
        ParametersPane1.add(labelPanel3);
        ParametersPane1.add(Box.createRigidArea(new Dimension(5,0)));
        ParametersPane1.add(TablePane1);
    
        ParametersPane2.add(labelPanel4);
        ParametersPane2.add(Box.createRigidArea(new Dimension(1,0)));
        ParametersPane2.add(TablePane2);

        ImageIcon upButtonIcon = createImageIcon("UP.gif");
        ImageIcon downButtonIcon = createImageIcon("down.gif");      
        
        JPanel panelUpDown = new JPanel();
        
        panelUpDown.setLayout(new BoxLayout(panelUpDown, BoxLayout.PAGE_AXIS));
        
        //JPanel panelSpecies12 = new JPanel();
        //panelSpecies12.setLayout(new GridLayout(1,2));
        String[] countOfMoleculesStrings = {"1","2","3","4","5","6","7","8"};
        SpinnerListModel countOfSpecies1SLM = new SpinnerListModel(countOfMoleculesStrings);
        SpinnerListModel countOfSpecies2SLM = new SpinnerListModel(countOfMoleculesStrings);
        SpinnerListModel countOfSpecies3SLM = new SpinnerListModel(countOfMoleculesStrings);
        
        JSpinner countSpecies1Spinner = new JSpinner(countOfSpecies1SLM);
        JSpinner countSpecies2Spinner = new JSpinner(countOfSpecies2SLM);
        JSpinner countSpecies3Spinner = new JSpinner(countOfSpecies3SLM);
        
        
        JPanel temp = new JPanel();
        //temp.setLayout(new GridLayout(5,1));
        temp.setLayout(new BoxLayout(temp, BoxLayout.PAGE_AXIS) );

        temp.add(countSpecies1Spinner);
        temp.add(Box.createRigidArea(new Dimension(0,10)));
        temp.add(countSpecies2Spinner);
        temp.add(Box.createRigidArea(new Dimension(0,10)));
        temp.add(countSpecies3Spinner);

        TablePane2.add(temp);
        
  //      panelUpDown.add(panelSpecies12);
  //      panelUpDown.add(panelSpecies3);
  //       ParametersPane2.add(panelUpDown);
 
    
        Border compound3;
        compound3 = BorderFactory.createCompoundBorder(
        		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
        ParametersPane.setBorder(compound3);
    
        JPanel DescriptionButtonPane = new JPanel();
        DescriptionButtonPane.setLayout(new FlowLayout());
    
        AddSameSpecies = new JButton("+ same species");
        AddAnotherSpecies = new JButton("+ diff species");
        Remove = new JButton("remove species");
    
        AddSameSpecies.setEnabled(false);
        AddAnotherSpecies.setEnabled(false);
        Remove.setEnabled(false);
    
        DescriptionButtonPane.add(AddSameSpecies);
        DescriptionButtonPane.add(AddAnotherSpecies);
        DescriptionButtonPane.add(Remove);
    
    
        JPanel DescriptionButtonPane2 = new JPanel();
        DescriptionButtonPane2.setLayout(new FlowLayout());
        
        Reset = new JButton("Reset all choices");
        AlterSimEnv = new JButton("Simulation Environment");
        Run = new JButton("Run Simulation");
    
    
        AlterSimEnv.setEnabled(false);
        Run.setEnabled(false);
        Reset.setEnabled(false);

        DescriptionButtonPane2.add(Reset);
        DescriptionButtonPane2.add(AlterSimEnv);
        DescriptionButtonPane2.add(Run);

        moleculeListDisplayTable.setIntercellSpacing(new Dimension(0,30));
        
        JComponent[] componentLeft = {AddSpecies,AddPotential,EditVariables,new JPanel()};
        JComponent[] componentRight = {listPane,list2Pane,ParametersPane,DescriptionButtonPane2};
        addComponentLeftRightRows(componentLeft,componentRight,gridbaglayout,panel1);
        panel1.setBorder(BorderFactory.createLineBorder(Color.black));
	
        Description = new JTextArea();
        Description.setEditable(false);
    
        this.add(new JScrollPane(Description), BorderLayout.CENTER);
  
	}

		
		protected static ImageIcon createImageIcon(String path) {
		    java.net.URL imgURL = ViewSpeciesSelection.class.getResource(path);
		    //error handling omitted for clarity...
		    return new ImageIcon(imgURL);
		}
		
	public JButton getAlterSimEnv() {
			return AlterSimEnv;
		}

		public void setAlterSimEnv(JButton alterSimEnv) {
			AlterSimEnv = alterSimEnv;
		}

	private void addComponentLeftRightRows(JComponent[] ComponentLeft,
        JComponent[] ComponentRight,
        GridBagLayout gridbag,
        Container container) {
	
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.NORTHEAST;
		int numLabels = ComponentLeft.length;

		for (int i = 0; i < numLabels; i++) {
			c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
			c.fill = GridBagConstraints.NONE;      //reset to default
			c.weightx = 0.0;   
			c.insets = new Insets(0,10,0,0);//reset to default
			container.add(ComponentLeft[i], c);

			c.gridwidth = GridBagConstraints.REMAINDER;     //end row
			c.fill = GridBagConstraints.HORIZONTAL;
			c.weightx = 0.0;
			c.insets = new Insets(0,10,3,0);
			container.add(ComponentRight[i], c);
		}
}
	
	public JButton getAddSpecies() {
		return AddSpecies;
	}



	public void setAddSpecies(JButton addSpecies) {
		AddSpecies = addSpecies;
	}



	public JButton getAddPotential() {
		return AddPotential;
	}



	public void setAddPotential(JButton addPotential) {
		AddPotential = addPotential;
	}



	public JButton getEditVariables() {
		return EditVariables;
	}



	public void setEditVariables(JButton editVariables) {
		EditVariables = editVariables;
	}



	public JButton getAddSameSpecies() {
		return AddSameSpecies;
	}



	public void setAddSameSpecies(JButton addSameSpecies) {
		AddSameSpecies = addSameSpecies;
	}



	public JButton getAddAnotherSpecies() {
		return AddAnotherSpecies;
	}



	public void setAddAnotherSpecies(JButton addAnotherSpecies) {
		AddAnotherSpecies = addAnotherSpecies;
	}



	public JButton getRemove() {
		return Remove;
	}



	public void setRemove(JButton remove) {
		Remove = remove;
	}



	public JButton getReset() {
		return Reset;
	}



	public void setReset(JButton reset) {
		Reset = reset;
	}



	public JButton getRun() {
		return Run;
	}



	public void setRun(JButton run) {
		Run = run;
	}

	
	public JTable getPotentialParamDispayTable() {
		return potentialParamDispayTable;
	}

	public void setPotentialParamDispayTable(JTable potentialParamDispayTable) {
		this.potentialParamDispayTable = potentialParamDispayTable;
	}

	public JTable getMoleculeListDisplayTable() {
		return moleculeListDisplayTable;
	}

	public void setMoleculeListDisplayTable(JTable moleculeListDisplayTable) {
		this.moleculeListDisplayTable = moleculeListDisplayTable;
	}

	public JTextArea getDescription() {
		return Description;
	}



	public void setDescription(JTextArea description) {
		Description = description;
	}



	public JList getSpeciesJList() {
		return speciesJList;
	}



	public void setSpeciesJList(JList speciesJList) {
		this.speciesJList = speciesJList;
	}



	public JList getPotentialJList() {
		return potentialJList;
	}



	public void setPotentialJList(JList potentialJList) {
		this.potentialJList = potentialJList;
	}

	public ListSelectionModel getSpeciesListDisplayLSM() {
		return speciesListDisplayLSM;
	}

	public void setSpeciesListDisplayLSM(ListSelectionModel speciesListDisplayLSM) {
		this.speciesListDisplayLSM = speciesListDisplayLSM;
	}

	public ListSelectionModel getPotentialListDisplayLSM() {
		return potentialListDisplayLSM;
	}

	public void setPotentialListDisplayLSM(
			ListSelectionModel potentialListDisplayLSM) {
		this.potentialListDisplayLSM = potentialListDisplayLSM;
	}

	public ModelPotParamDisplayTable getPotentialParamDisplayTM() {
		return potentialParamDisplayTM;
	}

	public void setPotentialParamDisplayTM(
			ModelPotParamDisplayTable potentialParamDisplayTM) {
		this.potentialParamDisplayTM = potentialParamDisplayTM;
	}

	public ModelMoleculeListDisplayTable_new getMoleculeListDisplayTM() {
		return moleculeListDisplayTM;
	}

	public void setMoleculeListDisplayTM(
			ModelMoleculeListDisplayTable_new moleculeListDisplayTM) {
		this.moleculeListDisplayTM = moleculeListDisplayTM;
	}

	public ListSelectionModel getPotentialParamDisplayLSM() {
		return potentialParamDisplayLSM;
	}

	public void setPotentialParamDisplayLSM(
			ListSelectionModel potentialParamDisplayLSM) {
		this.potentialParamDisplayLSM = potentialParamDisplayLSM;
	}

	public ListSelectionModel getMoleculeListDisplayLSM() {
		return moleculeListDisplayLSM;
	}

	public void setMoleculeListDisplayLSM(ListSelectionModel moleculeListDisplayLSM) {
		this.moleculeListDisplayLSM = moleculeListDisplayLSM;
	}
	
	
	
	public TableModelListener getPotentialParamDisplayTMListener() {
		return potentialParamDisplayTMListener;
	}

	public void setPotentialParamDisplayTMListener(
			TableModelListener potentialParamDisplayTMListener) {
		this.potentialParamDisplayTMListener = potentialParamDisplayTMListener;
	}

	public TableModelListener getMoleculeListDisplayTMListener() {
		return moleculeListDisplayTMListener;
	}

	public void setMoleculeListDisplayTMListener(
			TableModelListener moleculeListDisplayTMListener) {
		this.moleculeListDisplayTMListener = moleculeListDisplayTMListener;
	}

	public ListSelectionListener getSpeciesJListListener() {
		return speciesJListListener;
	}

	public void setSpeciesJListListener(ListSelectionListener speciesJListListener) {
		this.speciesJListListener = speciesJListListener;
	}

	public ListSelectionListener getPotentialJListListener() {
		return potentialJListListener;
	}

	public void setPotentialJListListener(
			ListSelectionListener potentialJListListener) {
		this.potentialJListListener = potentialJListListener;
	}

	
	
	public JButton getAddSpecies1() {
		return addSpecies1;
	}

	public void setAddSpecies1(JButton addSpecies1) {
		this.addSpecies1 = addSpecies1;
	}

	public JButton getAddSpecies2() {
		return addSpecies2;
	}

	public void setAddSpecies2(JButton addSpecies2) {
		this.addSpecies2 = addSpecies2;
	}

	public JButton getAddSpecies3() {
		return addSpecies3;
	}

	public void setAddSpecies3(JButton addSpecies3) {
		this.addSpecies3 = addSpecies3;
	}

	public JButton getRemoveSpecies1() {
		return removeSpecies1;
	}

	public void setRemoveSpecies1(JButton removeSpecies1) {
		this.removeSpecies1 = removeSpecies1;
	}

	public JButton getRemoveSpecies2() {
		return removeSpecies2;
	}

	public void setRemoveSpecies2(JButton removeSpecies2) {
		this.removeSpecies2 = removeSpecies2;
	}

	public JButton getRemoveSpecies3() {
		return removeSpecies3;
	}

	public void setRemoveSpecies3(JButton removeSpecies3) {
		this.removeSpecies3 = removeSpecies3;
	}

	public JTextField getSpeciesCount1() {
		return speciesCount1;
	}

	public void setSpeciesCount1(JTextField speciesCount1) {
		this.speciesCount1 = speciesCount1;
	}

	public JTextField getSpeciesCount2() {
		return speciesCount2;
	}

	public void setSpeciesCount2(JTextField speciesCount2) {
		this.speciesCount2 = speciesCount2;
	}

	public JTextField getSpeciesCount3() {
		return speciesCount3;
	}

	public void setSpeciesCount3(JTextField speciesCount3) {
		this.speciesCount3 = speciesCount3;
	}

	
	
	public void displayMoleculeList(){
		moleculeListDisplayTM.removeData();
		try{
		for(int i = 0;i < modelSpeciesSelection.getMoleculeListDisplayString().size() ; i++){
			moleculeListDisplayTM.setValueAt("Species " + Integer.toString(i + 1), i, 0);
			moleculeListDisplayTM.setValueAt(modelSpeciesSelection.getMoleculeListDisplayString().get(i), i, 1);
		}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	
	public void reset(){
		 
	}

	
	
	//Listener method to run the simulation when the Run button is pressed
	public void addSpeciesButtonListener(ActionListener mal) {
		AddSpecies.addActionListener(mal);
	}
	
	public void addSpeciesJListListener(ListSelectionListener mal) {
		setSpeciesJListListener(mal);
		speciesJList.addListSelectionListener(speciesJListListener);
	}
	
	public void removeSpeciesJListListener() {
		speciesJList.removeListSelectionListener(speciesJListListener);
	}
	
	
	public void addPotentialButtonListener(ActionListener mal) {
		AddPotential.addActionListener(mal);
	}
	
	public void addPotentialsJListListener(ListSelectionListener mal) {
		setPotentialJListListener(mal);
		potentialJList.addListSelectionListener(mal);
	}
	
	public void removePotentialJListListener() {
		potentialJList.removeListSelectionListener(potentialJListListener);
	}
	
	public void addEditVariablesButtonListener(ActionListener mal) {
		EditVariables.addActionListener(mal);
	}
	
	public void addPotParamDisplayTMListener(TableModelListener tal){
		setPotentialParamDisplayTMListener(tal);
		potentialParamDisplayTM.addTableModelListener(potentialParamDisplayTMListener);
	}
	
	public void removePotParamDisplayTMListener(){
		potentialParamDisplayTM.removeTableModelListener(potentialParamDisplayTMListener);
	}
	
	public void addMoleculeListDisplayTMListener(TableModelListener tal){
		setMoleculeListDisplayTMListener(tal);
		moleculeListDisplayTM.addTableModelListener(moleculeListDisplayTMListener);
		
	}
	
	public void removeMoleculeListDisplayTMListener(){
		moleculeListDisplayTM.removeTableModelListener(moleculeListDisplayTMListener);
	}
	
	public void addPotParamDisplayLSMListener(ListSelectionListener mal){
		potentialParamDisplayLSMListener = mal;
		potentialParamDisplayLSM.addListSelectionListener(potentialParamDisplayLSMListener);
	}
	
	public void removePotParamDisplayLSMListener(){
		potentialParamDisplayLSM.removeListSelectionListener(potentialParamDisplayLSMListener);
	}
	
	
	
	public void addMoleculeListDisplayLSMListener(ListSelectionListener mal){
		moleculeListDisplayLSMListener = mal;
		moleculeListDisplayLSM.addListSelectionListener(moleculeListDisplayLSMListener);
	}
	
	public void removeMoleculeListDisplayLSMListener(){
		moleculeListDisplayLSM.removeListSelectionListener(moleculeListDisplayLSMListener);
	}
	
	
	
	
	
	public void addAlterSimEnvButtonListener(ActionListener mal) {
		
		AlterSimEnv.addActionListener(mal);
	}
	
	public void addRunButtonListener(ActionListener mal) {
		Run.addActionListener(mal);
	}
	
	public void addRemoveButtonListener(ActionListener mal) {
		Remove.addActionListener(mal);
	}
	
	public void addSameSpeciesButtonListener(ActionListener mal) {
		AddSameSpecies.addActionListener(mal);
	}
	
	public void addAnotherSpeciesButtonListener(ActionListener mal) {
		AddAnotherSpecies.addActionListener(mal);
	}
	
	public void addResetButtonListener(ActionListener mal) {
		Reset.addActionListener(mal);
	}
	
	
	public void addAddSpecies1ButtonListener(ActionListener mal){
		addSpecies1.addActionListener(mal);
	}
	
	public void addAddSpecies2ButtonListener(ActionListener mal){
		addSpecies2.addActionListener(mal);
	}
	
	public void addAddSpecies3ButtonListener(ActionListener mal){
		addSpecies3.addActionListener(mal);
	}
	
	
	public void addRemoveSpecies1ButtonListener(ActionListener mal){
		removeSpecies1.addActionListener(mal);
	}
	
	public void addRemoveSpecies2ButtonListener(ActionListener mal){
		removeSpecies1.addActionListener(mal);
	}
	
	public void addRemoveSpecies3ButtonListener(ActionListener mal){
		removeSpecies1.addActionListener(mal);
	}
	
	
	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				
				ViewSpeciesSelection s = new ViewSpeciesSelection(new ModelSpeciesSelection(), new ModelSelectedSpecies());
				JFrame frame = new JFrame();
				frame.add(s);
				frame.setVisible(true);
				frame.setResizable(true);
				frame.setMinimumSize(new Dimension(830,800));
				frame.setBounds(0, 0,830,800);
				
			}
		});
    }


}
