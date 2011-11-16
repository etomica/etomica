package etomica.virial.GUI.containers;

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

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;

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






import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.virial.GUI.components.CreatePotentialCollections;
import etomica.virial.GUI.components.CreateSimulation;
import etomica.virial.GUI.components.PotentialCollectionFactory;
import etomica.virial.GUI.components.SimulationEnvironment;
import etomica.virial.GUI.components.SimulationEnvironmentObject;
import etomica.virial.GUI.models.CreateSpeciesDM_Alkane_SKS;
import etomica.virial.GUI.models.CreateSpeciesDM_Alkane_TRAPPE;
import etomica.virial.GUI.models.CreateSpeciesDM_CO2_2CLJQ;
import etomica.virial.GUI.models.CreateSpeciesDM_CO2_EMP2;
import etomica.virial.GUI.models.CreateSpeciesDM_CO2_Trappe;
import etomica.virial.GUI.models.CreateSpeciesDM_Ethane_2CLJQ;
import etomica.virial.GUI.models.CreateSpeciesDM_H2O_SPCE;
import etomica.virial.GUI.models.CreateSpeciesDM_LJ_2CLJQ;
import etomica.virial.GUI.models.CreateSpeciesDM_LJ_LJ;
import etomica.virial.GUI.models.CreateSpeciesDM_LJ_LJQ;
import etomica.virial.GUI.models.MixtureBuilderDM_ListingTable;
import etomica.virial.GUI.models.PotentialParamDM_Description;
import etomica.virial.GUI.models.CreateSpeciesDM_IFactory;


public class SpeciesSelectionUIView extends JPanel implements ActionListener, ListSelectionListener, TableModelListener{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//All JButtons which need to be accessed
	
	private JButton AddSpecies;
	private JButton AddPotential;
	private JButton EditVariables;
	
	private JButton AddSameSpecies;
	private JButton AddAnotherSpecies;
	private JButton Remove;
	
	private JButton Reset;
	private JButton OtherParamValues;
	private JButton Run;
	
	//This needs to be removed once this code is ready to be added to the main code
	
	
	
	
	private String[] IntialList = {"---No species selected---","Press \'Add\' to select a species"};
	private String[] IntialPotentialList = {"---No potentials selected---","Press \'Add\' to select a potential"};
	
	private String[] SpeciesList = {"LJ","CO2","Methanol","Ethanol","Methane","Ethane","Propane","Higher n-Alkanes","Naphthalene","Water"};
	
	@SuppressWarnings("rawtypes")
	private Class[] LJPotentialClassList = {CreateSpeciesDM_LJ_LJ.class,CreateSpeciesDM_LJ_LJQ.class,CreateSpeciesDM_LJ_2CLJQ.class};
	
	@SuppressWarnings("rawtypes")
	private Class[] CO2PotentialClassList = {CreateSpeciesDM_CO2_2CLJQ.class,CreateSpeciesDM_CO2_EMP2.class,CreateSpeciesDM_CO2_Trappe.class};
	
	@SuppressWarnings("rawtypes")
	private Class[] MethanePotentialClassList = {CreateSpeciesDM_Alkane_TRAPPE.class};
	
	@SuppressWarnings("rawtypes")
	private Class[] EthanePotentialClassList = {CreateSpeciesDM_Alkane_TRAPPE.class, CreateSpeciesDM_Alkane_SKS.class, CreateSpeciesDM_Ethane_2CLJQ.class};
	
	@SuppressWarnings("rawtypes")
	private Class[] PropanePotentialClassList = {CreateSpeciesDM_Alkane_TRAPPE.class,CreateSpeciesDM_Alkane_SKS.class};
	
	@SuppressWarnings("rawtypes")
	private Class[] nAlkanePotentialClassList = {CreateSpeciesDM_Alkane_TRAPPE.class, CreateSpeciesDM_Alkane_SKS.class};
	
	@SuppressWarnings("rawtypes")
	private Class[] WaterPotentialClassList = {CreateSpeciesDM_H2O_SPCE.class};
	
	
	
	//For each species...
	private String PotentialType = null;
	private int NoOfVariables;
	private int[] SpeciesIndex = new int[2];
	private int[] PotentialsIndex = new int[2];
	
	private int Alkane1Spheres = 4;
	private int Alkane2Spheres = 4;
	private AlkaneSpheresUIView NSpheres1;
	private AlkaneSpheresUIView NSpheres2;
	private JFrame NAlkaneFrame1;
	private JFrame NAlkaneFrame2;
	
	//Can be separated From Species Class!!
	private SimulationEnvironmentUIView SimulationEnvParam;
	private JFrame OtherParamViewFrame;
	private SimulationEnvironmentObject OtherParamObject;
	private SimulationEnvironment SimENV;
	private PotentialCollectionFactory PObject;
	private CreateSimulation simulation;
	private Console console;
	private AlertMessageUIView MessageAlert; 
	private JFrame MessageFrame;
	



	//Parameter Tables
	private JTable table1;
	private JTable table2;

	//Description
	private JTextArea Description; 
	
	//JList 
	private JList SpeciesJList;
	private JList PotentialJList;
	
	private MixtureBuilderDM_ListingTable LJP1;
	private MixtureBuilderDM_ListingTable LJP2;
	
	private ListSelectionModel cellSelectionModel;
	private ListSelectionModel cellSelectionModel2;
	private ListSelectionModel SpeciesListSelectionModel;
	private ListSelectionModel PotentialListSelectionModel;
	
	private CreateSpeciesDM_IFactory potential1;
	private CreateSpeciesDM_IFactory potential2;
	private int TotalMoleculeCount = 0;
	private int Molecule1Count = 0;
	private int Molecule2Count = 0;
	
	private boolean AddFlag = false;
	private boolean ChangedSpecies = false;
	
	//private SpeciesList sList;
	private ArrayList<String> sList;
	
	public SpeciesSelectionUIView(){
		super();
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
		//AddSpecies.addActionListener(this);
	
		SpeciesJList = new JList(IntialList);
		SpeciesJList.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		SpeciesListSelectionModel = SpeciesJList.getSelectionModel();
		SpeciesJList.setLayoutOrientation(JList.HORIZONTAL_WRAP);
		SpeciesJList.setVisibleRowCount(-1);
		SpeciesJList.addListSelectionListener(this);
		SpeciesJList.setBorder(BorderFactory.createLineBorder(Color.black));
	
 
		JScrollPane SpeciesListScroller = new JScrollPane(SpeciesJList);
		SpeciesListScroller.setPreferredSize(new Dimension(400,100));
		SpeciesListScroller.setBorder(BorderFactory.createLoweredBevelBorder());
    
		JPanel listPane = new JPanel();				
		listPane.setLayout(new BoxLayout(listPane, BoxLayout.PAGE_AXIS));
    
		JLabel label = new JLabel("Species");
		label.setForeground(Color.WHITE);
		JPanel labelPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		labelPanel.setBackground(new Color(0,78,152));
		labelPanel.setBounds(listPane.getBounds().x, listPane.getBounds().y, 400, label.getHeight());
   
		//labelPanel.setBackground(Color.green);
		label.setLabelFor(SpeciesJList);
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
		AddPotential.addActionListener(this);
   
		PotentialJList = new JList(IntialPotentialList);
	
		PotentialJList.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		PotentialJList.setLayoutOrientation(JList.HORIZONTAL_WRAP);
		PotentialJList.setVisibleRowCount(-1);
		PotentialJList.setEnabled(false);
		PotentialJList.addListSelectionListener(this);
 
		JScrollPane PotentialListScroller = new JScrollPane(PotentialJList);
		PotentialListScroller.setPreferredSize(new Dimension(400,100));
		PotentialListScroller.setBorder(BorderFactory.createLoweredBevelBorder());
  
		JPanel list2Pane = new JPanel();
		list2Pane.setLayout(new BoxLayout(list2Pane, BoxLayout.PAGE_AXIS));
   
		JLabel label2 = new JLabel("Potentials");
		label2.setForeground(Color.WHITE);
		label2.setLabelFor(PotentialJList);
    
		JPanel labelPanel2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		labelPanel2.setBackground(new Color(0,78,152));
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
		EditVariables.addActionListener(this);
    
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
    
		LJP1 = new MixtureBuilderDM_ListingTable();
		LJP2 = new MixtureBuilderDM_ListingTable();
    
		
		table1 = new JTable();
		table1.setModel(LJP1);
		table1.setEnabled(false);
   
		cellSelectionModel = table1.getSelectionModel();
		cellSelectionModel.addListSelectionListener(this);
		
		TableColumn column1 = null;
		for (int i = 0; i < 2; i++) {
			column1 = table1.getColumnModel().getColumn(i);
			column1.setPreferredWidth(150);
		}
    
        table2 = new JTable();
        table2.setModel(LJP2);
        table2.setEnabled(false);
    

		cellSelectionModel2 = table2.getSelectionModel();
		cellSelectionModel2.addListSelectionListener(this);
    
        
        TableColumn column2 = null;
        for (int i = 0; i < 2; i++) {
        	column2 = table2.getColumnModel().getColumn(i);
        	column2.setPreferredWidth(150);
       
        }
    
        TablePane1.add(table1);
        TablePane2.add(table2);
    
    
        JLabel label3 = new JLabel("Parameters");
        label3.setForeground(Color.WHITE);
        JPanel labelPanel3 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        labelPanel3.setBackground(new Color(0,78,152));
    
        /*labelPanel3.setBounds(ParametersPane.getBounds().x, ParametersPane.getBounds().y, 400, label3.getHeight());*/
        labelPanel3.setBounds(ParametersPane1.getBounds().x, ParametersPane2.getBounds().y, 400, label3.getHeight());
        labelPanel3.add(label3);
    
        JLabel label4 = new JLabel("Molecule List");
        label4.setForeground(Color.WHITE);
        JPanel labelPanel4 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        labelPanel4.setBackground(new Color(0,78,152));
        labelPanel4.setBounds(ParametersPane.getBounds().x + labelPanel3.getWidth(), ParametersPane.getBounds().y, 400, label4.getHeight());
        labelPanel4.add(label4);
    
        JPanel labelPane = new JPanel();
        labelPane.add(labelPanel3);
        labelPane.add(labelPanel4);
   
    
        ParametersPane1.add(labelPanel3);
        ParametersPane1.add(Box.createRigidArea(new Dimension(5,0)));
        ParametersPane1.add(TablePane1);
    
        ParametersPane2.add(labelPanel4);
        ParametersPane2.add(Box.createRigidArea(new Dimension(5,0)));
        ParametersPane2.add(TablePane2);
    
        Border compound3;
        compound3 = BorderFactory.createCompoundBorder(
        		BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder());
        ParametersPane.setBorder(compound3);
    
        JPanel DescriptionButtonPane = new JPanel();
        DescriptionButtonPane.setLayout(new FlowLayout());
    
        AddSameSpecies = new JButton("+ same species");
        AddAnotherSpecies = new JButton("+ diff species");
        Remove = new JButton("remove species");
    
        Remove.addActionListener(this);
        AddSameSpecies.addActionListener(this);
        AddAnotherSpecies.addActionListener(this);
    
        DescriptionButtonPane.add(AddSameSpecies);
        DescriptionButtonPane.add(AddAnotherSpecies);
        DescriptionButtonPane.add(Remove);
    
    
        JPanel DescriptionButtonPane2 = new JPanel();
        DescriptionButtonPane2.setLayout(new FlowLayout());
        
        Reset = new JButton("Reset all choices");
        OtherParamValues = new JButton("Simulation Environment");
        Run = new JButton("Run Simulation");
    
    
        Reset.addActionListener(this);
        OtherParamValues.addActionListener(this);
        Run.addActionListener(this);
        DescriptionButtonPane2.add(Reset);
        DescriptionButtonPane2.add(OtherParamValues);
        DescriptionButtonPane2.add(Run);

        JComponent[] componentLeft = {AddSpecies,AddPotential,EditVariables,new JPanel(),new JPanel()};
        JComponent[] componentRight = {listPane,list2Pane,ParametersPane,DescriptionButtonPane,DescriptionButtonPane2};
        addComponentLeftRightRows(componentLeft,componentRight,gridbaglayout,panel1);
        panel1.setBorder(BorderFactory.createLineBorder(Color.black));
	
        Description = new JTextArea();
        Description.setEditable(false);
    
        this.add(new JScrollPane(Description), BorderLayout.CENTER);
  
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
	
	public void resetToAddAnotherSpecies(){
		SpeciesJList.removeListSelectionListener(this);
		SpeciesJList.clearSelection();
		SpeciesJList.setListData(IntialList);
		
		cellSelectionModel.removeListSelectionListener(this);
		cellSelectionModel.clearSelection();
		
		PotentialType = null;
		PotentialJList.removeListSelectionListener(this);
		PotentialJList.clearSelection();
		PotentialJList.setListData(IntialPotentialList);
		if(!Remove.isEnabled()){
			Remove.setEnabled(true);
		}
		if(!AddSameSpecies.isEnabled()){
			AddSameSpecies.setEnabled(true);
		}
		if(!AddAnotherSpecies.isEnabled()){
			AddAnotherSpecies.setEnabled(true);
		}
		LJP1.removeTableModelListener(this);
		LJP1.removeData();
		SpeciesJList.addListSelectionListener(this);
		PotentialJList.addListSelectionListener(this);
		cellSelectionModel.addListSelectionListener(this);
		AddPotential.setEnabled(false);
		EditVariables.setEnabled(false);
		
		
		displayMoleculeList();
	}
	
	public void resetRemovingPotential2(){
		SpeciesJList.removeListSelectionListener(this);
		SpeciesJList.clearSelection();
		SpeciesJList.setListData(SpeciesList);
		SpeciesJList.setSelectedIndex(SpeciesIndex[0]);
		
		cellSelectionModel.removeListSelectionListener(this);
		cellSelectionModel.clearSelection();
		PotentialJList.removeListSelectionListener(this);
		PotentialJList.clearSelection();
		switch(SpeciesIndex[0]){
			case 0:
				PotentialType = "LJ";
				setPotentialList(LJPotentialClassList,PotentialType);
				break;
			case 1:
				PotentialType = "CO2";
				setPotentialList(CO2PotentialClassList,PotentialType);
				break;
			case 2:
				PotentialType = "Methanol";
				break;
			case 3:
				PotentialType = "Ethanol";
				break;
			case 4:
				PotentialType = "Methane";
				setPotentialList(MethanePotentialClassList,PotentialType);
				break;
				
			case 5:
				PotentialType = "Ethane";
				setPotentialList(EthanePotentialClassList,PotentialType);
				break;
			case 6:
				PotentialType = "Propane";
				setPotentialList(PropanePotentialClassList,PotentialType);
				break;
			case 7:
				PotentialType = "Napthalene";
				break;
			case 9:
				PotentialType = "Water";
				setPotentialList(WaterPotentialClassList,PotentialType);
				break;
		}
		
		//PotentialJList.setListData(IntialPotentialList);
		PotentialJList.setSelectedIndex(PotentialsIndex[0]);
		if(!Remove.isEnabled()){
			Remove.setEnabled(true);
		}
		if(!AddSameSpecies.isEnabled()){
			AddSameSpecies.setEnabled(true);
		}
		if(!AddAnotherSpecies.isEnabled()){
			AddAnotherSpecies.setEnabled(true);
		}
		LJP1.removeTableModelListener(this);
		LJP1.removeData();
		String[][] ParamAndValue = potential1.getParamAndValues();
		SpeciesJList.addListSelectionListener(this);
		PotentialJList.addListSelectionListener(this);
		cellSelectionModel.addListSelectionListener(this);
		AddPotential.setEnabled(false);
		EditVariables.setEnabled(false);
		
		
		for(int i = 0;i <ParamAndValue.length;i++){
			LJP1.setValueAt(ParamAndValue[i][0],i,0);
			LJP1.setValueAt(ParamAndValue[i][1],i,1);
		}
		
		displayMoleculeList();
		AddFlag = false;
	}
	
	public void displayMoleculeList(){
		LJP2.removeData();
		ArrayList<String> DisplayContents = sList;
		for(int i = 0;i < DisplayContents.size() ; i++){
			LJP2.setValueAt("Molecule " + Integer.toString(i + 1), i, 0);
			LJP2.setValueAt(DisplayContents.get(i), i, 1);
		}
	}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void setPotentialList(Class[] potentialList,String potentialType){
		String[] s = new String[potentialList.length];
        for (int i=0; i<s.length; i++) {
		   try {
			try {
				s[i] = (String)potentialList[i].getMethod("getCustomName", new Class[0]).invoke(potentialList[i].getConstructor().newInstance(new Object[0]),new Object[0]);
			} catch (InstantiationException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		} catch (IllegalArgumentException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (SecurityException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalAccessException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (InvocationTargetException e1) {
			// TODO Auto-generated catch blockPotentialJList
			e1.printStackTrace();
		} catch (NoSuchMethodException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		}
        try{
        	PotentialType = potentialType;
        	PotentialJList.setListData(s);
        	
        }
        catch (Exception E){
			E.printStackTrace();
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



	public String[] getIntialList() {
		return IntialList;
	}



	public void setIntialList(String[] intialList) {
		IntialList = intialList;
	}



	public String[] getIntialPotentialList() {
		return IntialPotentialList;
	}



	public void setIntialPotentialList(String[] intialPotentialList) {
		IntialPotentialList = intialPotentialList;
	}



	public String[] getSpeciesList() {
		return SpeciesList;
	}



	public void setSpeciesList(String[] speciesList) {
		SpeciesList = speciesList;
	}



	public Class[] getLJPotentialClassList() {
		return LJPotentialClassList;
	}



	public void setLJPotentialClassList(Class[] lJPotentialClassList) {
		LJPotentialClassList = lJPotentialClassList;
	}



	public Class[] getCO2PotentialClassList() {
		return CO2PotentialClassList;
	}



	public void setCO2PotentialClassList(Class[] cO2PotentialClassList) {
		CO2PotentialClassList = cO2PotentialClassList;
	}



	public Class[] getMethanePotentialClassList() {
		return MethanePotentialClassList;
	}



	public void setMethanePotentialClassList(Class[] methanePotentialClassList) {
		MethanePotentialClassList = methanePotentialClassList;
	}



	public Class[] getEthanePotentialClassList() {
		return EthanePotentialClassList;
	}



	public void setEthanePotentialClassList(Class[] ethanePotentialClassList) {
		EthanePotentialClassList = ethanePotentialClassList;
	}



	public Class[] getPropanePotentialClassList() {
		return PropanePotentialClassList;
	}



	public void setPropanePotentialClassList(Class[] propanePotentialClassList) {
		PropanePotentialClassList = propanePotentialClassList;
	}



	public String getPotentialType() {
		return PotentialType;
	}



	public void setPotentialType(String potentialType) {
		PotentialType = potentialType;
	}



	public int getNoOfVariables() {
		return NoOfVariables;
	}



	public void setNoOfVariables(int noOfVariables) {
		NoOfVariables = noOfVariables;
	}



	public JTable getTable1() {
		return table1;
	}



	public void setTable1(JTable table1) {
		this.table1 = table1;
	}



	public JTable getTable2() {
		return table2;
	}



	public void setTable2(JTable table2) {
		this.table2 = table2;
	}



	public JTextArea getDescription() {
		return Description;
	}



	public void setDescription(JTextArea description) {
		Description = description;
	}



	public JList getSpeciesJList() {
		return SpeciesJList;
	}



	public void setSpeciesJList(JList speciesJList) {
		SpeciesJList = speciesJList;
	}



	public JList getPotentialJList() {
		return PotentialJList;
	}



	public void setPotentialJList(JList potentialJList) {
		PotentialJList = potentialJList;
	}



	public MixtureBuilderDM_ListingTable getLJP1() {
		return LJP1;
	}



	public void setLJP1(MixtureBuilderDM_ListingTable lJP1) {
		LJP1 = lJP1;
	}



	public MixtureBuilderDM_ListingTable getLJP2() {
		return LJP2;
	}



	public void setLJP2(MixtureBuilderDM_ListingTable lJP2) {
		LJP2 = lJP2;
	}



	public ListSelectionModel getCellSelectionModel() {
		return cellSelectionModel;
	}



	public void setCellSelectionModel(ListSelectionModel cellSelectionModel) {
		this.cellSelectionModel = cellSelectionModel;
	}



	public CreateSpeciesDM_IFactory getPotential1() {
		return potential1;
	}



	public void setPotential1(CreateSpeciesDM_IFactory potential1) {
		this.potential1 = potential1;
	}



	public CreateSpeciesDM_IFactory getPotential2() {
		return potential2;
	}



	public void setPotential2(CreateSpeciesDM_IFactory potential2) {
		this.potential2 = potential2;
	}



	public boolean isAddFlag() {
		return AddFlag;
	}



	public void setAddFlag(boolean addFlag) {
		AddFlag = addFlag;
	}



	public boolean isChangedSpecies() {
		return ChangedSpecies;
	}



	public void setChangedSpecies(boolean changedSpecies) {
		ChangedSpecies = changedSpecies;
	}


/*
	public SpeciesList getsList() {
		return sList;
	}



	public void setsList(SpeciesList sList) {
		this.sList = sList;
	}
*/
	public ArrayList<String> getsList() {
		return sList;
	}



	public void setsList(ArrayList<String> sList) {
		this.sList = sList;
	}
	
	
public void actionPerformed(ActionEvent e){
		
		if(e.getSource().equals(AddSpecies)){
			
			
				SpeciesJList.setListData(SpeciesList);
				PotentialJList.setEnabled(true);
				AddPotential.setEnabled(true);
				
				//This AddFlag is true when a Different Species is to be added to the list
				//When AddFlag is false, we re changing the species selection for potential 1, thats why we give a new instance
				//of the the sList ArrayList
				if(!AddFlag){
					//sList = new SpeciesList();
					sList = new ArrayList<String>();
				}
				
				
			
		}
		
		if(e.getSource().equals(AddPotential)){
			
				ChangedSpecies = true;
				
				//This will store the values of the indexes selected
				if(AddFlag){
					SpeciesIndex[1] = SpeciesJList.getSelectedIndex();
				}
				else{
					SpeciesIndex[0] = SpeciesJList.getSelectedIndex();
				}
				
				if(SpeciesJList.getSelectedIndex() == 0){
					Description.setText(" ");
					Description.append("We now create a single LJ Molecule\n");
					setPotentialList(LJPotentialClassList,"LJ");
				}
				if(SpeciesJList.getSelectedIndex() == 1){
					Description.setText(" ");
					Description.append("We now create a single CO2 Molecule\n");
					setPotentialList(CO2PotentialClassList,"CO2");
				}
				if(SpeciesJList.getSelectedIndex() == 4){
					Description.setText(" ");
					Description.append("We now create a single methane Molecule\n");
					setPotentialList(MethanePotentialClassList,"Methane");

				}
			
				if(SpeciesJList.getSelectedIndex() == 5){
					Description.setText(" ");
					Description.append("We now create a single ethane Molecule\n");
					setPotentialList(EthanePotentialClassList,"Ethane");

				}
			
				if(SpeciesJList.getSelectedIndex() == 6){
					Description.setText(" ");
					Description.append("We now create a single propane Molecule\n");
					setPotentialList(PropanePotentialClassList,"Propane");

				}
				
				if(SpeciesJList.getSelectedIndex() == 7){
					Description.setText(" ");
					Description.append("We now create a single alkane Molecule\n");
					setPotentialList(nAlkanePotentialClassList,"n-Alkane");

				}
				
				if(SpeciesJList.getSelectedIndex() == 9){
					Description.setText(" ");
					Description.append("We now create a single Water Molecule\n");
					setPotentialList(WaterPotentialClassList,"Water");

				}
					}
		
		 
		if(e.getSource().equals(OtherParamValues)){
			OtherParamViewFrame = new JFrame("Simulation Environment");
			SimulationEnvironment SimENV = SimulationEnvironment.getInstance();
			if(OtherParamObject == null){
				OtherParamObject = new SimulationEnvironmentObject(SimENV.getTemperature(),SimENV.noOfSteps, potential1,potential2);
			}
			if(SimulationEnvParam == null){
			SimulationEnvParam = new SimulationEnvironmentUIView(OtherParamObject);}
			OtherParamViewFrame.add(SimulationEnvParam);
			SimulationEnvParam.getCloseWindow().addActionListener(this);
			SimulationEnvParam.getSaveValues().addActionListener(this);
			OtherParamViewFrame.setMinimumSize(new Dimension(500,250));
			OtherParamViewFrame.setBounds(this.getParent().getParent().getWidth(), 0, 0, 0);
			
			OtherParamViewFrame.setVisible(true);
			if(NAlkaneFrame1 != null){
				NAlkaneFrame1.setVisible(true);
			}
			if(NAlkaneFrame2 != null){
				NAlkaneFrame2.setVisible(true);
			}
		}
		
		if(SimulationEnvParam!= null){
			if(e.getSource().equals(SimulationEnvParam.getCloseWindow())){
				OtherParamViewFrame.setVisible(false);
			}
			
			if(e.getSource().equals(SimulationEnvParam.getSaveValues())){
				
				if(SimENV == null){
					SimENV = SimulationEnvironment.getInstance();
				}
				SimENV.setTemperature(Double.parseDouble(SimulationEnvParam.getTemperatureField().getText()));
				SimENV.setNoOfSteps(Integer.parseInt(SimulationEnvParam.getNoOfStepsField().getText()));
				if(SimulationEnvParam.getSigmaHSRefFieldA().getText() != "0.0"){
					potential1.setParameter("SIGMAHSREF",SimulationEnvParam.getSigmaHSRefFieldA().getText());}
				if(potential2 != null){
					if(SimulationEnvParam.getSigmaHSRefFieldB().getText() != "0.0"){
						potential2.setParameter("SIGMAHSREF",SimulationEnvParam.getSigmaHSRefFieldB().getText());
					}
				}
				OtherParamObject.setTemperature(Double.parseDouble(SimulationEnvParam.getTemperatureField().getText()));
				OtherParamObject.setNoOfSteps(Integer.parseInt(SimulationEnvParam.getNoOfStepsField().getText()));
				OtherParamObject.setSigmaHSRefA(Double.parseDouble(SimulationEnvParam.getSigmaHSRefFieldA().getText()));
				OtherParamObject.setSigmaHSRefB(Double.parseDouble(SimulationEnvParam.getSigmaHSRefFieldB().getText()));
				OtherParamViewFrame.setVisible(false);
			}
			
		}
		
		if(MessageAlert != null){
			if(e.getSource().equals(MessageAlert.getCloseWindow())){
				MessageFrame.setVisible(false);
				//resetToAddAnotherSpecies();
				resetRemovingPotential2();
				
			}
		}
		
		if(NSpheres1!= null){
			
			if(e.getSource().equals(NSpheres1.getCloseWindow())){
				Alkane1Spheres = Integer.parseInt(NSpheres1.getNoOfSpheres().getText());
				potential1.setParameter("NUMBER", Integer.toString(Alkane1Spheres));
				if(OtherParamObject != null){
					OtherParamObject.setAlkane1Spheres(Alkane1Spheres);
				}
				NAlkaneFrame1.setVisible(false);
			}
			
			if(e.getSource().equals(NSpheres1.getSaveValues())){
				Alkane1Spheres = Integer.parseInt(NSpheres1.getNoOfSpheres().getText());
				potential1.setParameter("NUMBER", Integer.toString(Alkane1Spheres));
				if(OtherParamObject != null){
					
					OtherParamObject.setAlkane1Spheres(Alkane1Spheres);
					
				}
				NAlkaneFrame1.setVisible(false);
			}
			
		}
		
		if(NSpheres2!= null){
			if(e.getSource().equals(NSpheres2.getCloseWindow())){
				Alkane2Spheres = Integer.parseInt(NSpheres2.getNoOfSpheres().getText());
				potential2.setParameter("NUMBER", Integer.toString(Alkane2Spheres));
				if(OtherParamObject != null){
					
					OtherParamObject.setAlkane2Spheres(Alkane2Spheres);
					
				}
				NAlkaneFrame2.setVisible(false);
			}
			
			if(e.getSource().equals(NSpheres2.getSaveValues())){
				Alkane2Spheres = Integer.parseInt(NSpheres2.getNoOfSpheres().getText());
				potential2.setParameter("NUMBER", Integer.toString(Alkane2Spheres));
				if(OtherParamObject != null){
					
					OtherParamObject.setAlkane2Spheres(Alkane2Spheres);
					
				}
				NAlkaneFrame2.setVisible(false);
			}
			
		}
		
		if(e.getSource().equals(Run)){

			simulation = new CreateSimulation(potential1, potential2);
			if(SimENV == null){
				SimENV = SimulationEnvironment.getInstance();
			}
			
			/*if(console==null){
				console = new Console();}*/
	
			try {
				if(OtherParamObject == null){
					OtherParamObject = new SimulationEnvironmentObject(SimENV.getTemperature(),SimENV.getNoOfSteps(),potential1,potential2);
				}	
				OtherParamObject.setSpeciesA(Molecule1Count);
				OtherParamObject.setSpeciesB(Molecule2Count);
				OtherParamObject.setnPoints();
				OtherParamObject.calculateSystemSigmaHSRef(potential1, potential2);
				
				
				CreatePotentialCollections checkSpecies = new CreatePotentialCollections();
				PObject = checkSpecies.checkIfCompatible(potential1,potential2,OtherParamObject);
				if(PObject == null){
					if(potential2 == null){
						simulation.runSimulation(OtherParamObject);
					}
				}else{
					simulation.runSimulation(OtherParamObject,PObject);
				}
				
			} catch (NoSuchMethodException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		
		
		if(e.getSource().equals(EditVariables)){
			//table1.setEnabled(true);
			LJP1.addTableModelListener(this);
			table2.setEnabled(true);
			LJP2.addTableModelListener(this);	
		}
		
		if(e.getSource().equals(Reset)){
			SpeciesJList.removeListSelectionListener(this);
			SpeciesJList.clearSelection();
			SpeciesJList.setListData(IntialList);
			
			cellSelectionModel.removeListSelectionListener(this);
			cellSelectionModel.clearSelection();
			
			PotentialType = null;
			PotentialJList.removeListSelectionListener(this);
			PotentialJList.clearSelection();
			PotentialJList.setListData(IntialPotentialList);
			if(!Remove.isEnabled()){
				Remove.setEnabled(true);
			}
			if(!AddSameSpecies.isEnabled()){
				AddSameSpecies.setEnabled(true);
			}
			if(!AddAnotherSpecies.isEnabled()){
				AddAnotherSpecies.setEnabled(true);
			}
			LJP1.removeTableModelListener(this);
			LJP1.removeData();
			LJP2.removeTableModelListener(this);
			LJP2.removeData();
			SpeciesJList.addListSelectionListener(this);
			PotentialJList.addListSelectionListener(this);
			cellSelectionModel.addListSelectionListener(this);
			AddFlag = false;
			
		}
		
		if(e.getSource().equals(Remove)){
			//if(sList.getId() > 1){
			if(sList.size() > 1){
				//LJP2.removeTableModelListener(this);
				boolean Potential2ExistsFlag = false;
				String[][] ParameterArray;
				if(!AddSameSpecies.isEnabled()){
					AddSameSpecies.setEnabled(true);
				}
				if(!AddAnotherSpecies.isEnabled()){
					AddAnotherSpecies.setEnabled(true);
				}
				//sList.removeSpecies();
				if(potential2 != null){
					if(sList.get(TotalMoleculeCount - 1) == potential2.getMoleculeDisplayName()){
						for(int i = 0;i<TotalMoleculeCount - 1;i++){
							if(sList.get(i) == potential2.getMoleculeDisplayName()){
								sList.remove(TotalMoleculeCount - 1);
								TotalMoleculeCount--;
								Molecule2Count--;
								displayMoleculeList();
								Potential2ExistsFlag = true;
								break;
							}
						
						}
						if(!Potential2ExistsFlag){
							sList.remove(TotalMoleculeCount - 1);
							TotalMoleculeCount--;
							Molecule2Count--;
							resetRemovingPotential2();
						
							//displayMoleculeList();
						}
					}
					else{
						sList.remove(TotalMoleculeCount - 1);
						TotalMoleculeCount--;
						Molecule1Count--;
						displayMoleculeList();
					}
				}
				else{
					sList.remove(TotalMoleculeCount - 1);
					TotalMoleculeCount--;
					Molecule1Count--;
					displayMoleculeList();
				}
				
				//TotalMoleculeCount--;
				//displayMoleculeList();
				Description.append("You ve just removed the last added species!!\n");
				LJP2.addTableModelListener(this);
				
			}
			//if(sList.getId() <= 1){
			if(sList.size() <= 1){
				//sList.removeSpecies();
					//potential2 = null;
					Description.append("You cannot remove the remaining molecule!!\n");
					Remove.setEnabled(false);
					AddFlag = false;
					
					//AddSameSpecies.setEnabled(false);
			}
		}
		if(e.getSource().equals(AddSameSpecies)){
			//if(sList.getId() < 9){
			if(sList.size() < 9){
				if(!Remove.isEnabled()){
					Remove.setEnabled(true);
				}
				//ParameterMapping clone;
				//clone = (ParameterMapping)potential1.clone();
				
				//sList.addSpecies(clone);
				sList.add(potential1.getMoleculeDisplayName());
				TotalMoleculeCount++;
				Molecule1Count++;
				Description.append("1 more Molecule added!\n");
				//Description.append("So far " + Integer.toString(sList.getId())+ " molecules in the sim box\n");
				Description.append("So far " + Integer.toString(TotalMoleculeCount)+ " molecules in the sim box\n");
				displayMoleculeList();
				//if(sList.getId() == 8){
				if(sList.size() == 8){
					Description.append("Sorry! You cannot add any more molecules!!\n");
					AddSameSpecies.setEnabled(false);
					AddAnotherSpecies.setEnabled(false);
				}
			}
			
		}
		
		if(e.getSource().equals(AddAnotherSpecies)){
			if(!AddFlag){
				resetToAddAnotherSpecies();
				
				AddFlag = true;
				if(!Remove.isEnabled()){
					Remove.setEnabled(true);
				}
			}
			else{
				
				//ParameterMapping clone = (ParameterMapping)potential2.clone();
				//sList.addSpecies(clone);
				sList.add(potential2.getMoleculeDisplayName());
				TotalMoleculeCount++;
				Molecule2Count++;
				Description.append("So far " + Integer.toString(TotalMoleculeCount)+ " molecules in the sim box\n");
				Description.append("1 more Molecule added!\n");
				//Description.append("So far " + Integer.toString(sList.getId())+ " molecules in the sim box\n");
				
				displayMoleculeList();
				//if(sList.getId() == 8){
				if(sList.size() == 8){
					Description.append("Sorry! You cannot add any more molecules!!\n");
					AddAnotherSpecies.setEnabled(false);
					AddSameSpecies.setEnabled(false);
					
				}
			
			
			}
			
		}
	}
	
	
	
	@SuppressWarnings("unchecked")
	public void valueChanged(ListSelectionEvent e){
		if(e.getSource().equals(SpeciesJList)){
			
			if(!e.getValueIsAdjusting()){
				/*for(int i = 0;i < 8;i++){
					if(sList.getObject(i) != null && AddFlag == false){
						sList.removeSpeciesAtIndex(i);
					}
				}*/
				if(sList.size() > 1 && AddFlag == false){
					sList.clear();
					LJP2.removeData();
					LJP1.removeData();
				}
				ChangedSpecies = false;
				
				
				PotentialType = null;
				PotentialJList.setListData(IntialPotentialList);
				
				
				
			}
			
		}
		
		
		
			if(e.getSource().equals(PotentialJList)){
				if(!e.getValueIsAdjusting()){
					
					JList list = (JList)e.getSource();
					if(ChangedSpecies){
					//Declarations within the scope of this 
					if(AddFlag){
						PotentialsIndex[1] = list.getSelectedIndex();}
					else{
						PotentialsIndex[0] = list.getSelectedIndex();
					}
					Class[] potentialList = null;
					String[] ParameterArray = null;
					
					String[] PotentialSite = null;
					String[][] ParamAndValue = null;
					Object obj[] = new Object[1];
					
					
					if(PotentialType == "LJ"){
						potentialList = LJPotentialClassList;
						table1.setEnabled(true);
						table2.setEnabled(true);
						if(list.getSelectedIndex() == 0){
							Description.setText(" ");
							Description.append("We now create a single LJ Molecule\n");
							Description.append("With a Spherical-2-Body Potential\n");
							EditVariables.setEnabled(true);
							
							}
						
						if(list.getSelectedIndex() == 1){
							Description.setText(" ");
							Description.append("We now create a single LJ Molecule\n");
							Description.append("With a Spherical-2-Body-With-Quad Potential\n");
							EditVariables.setEnabled(true);
						}
						if(list.getSelectedIndex() == 2){
							Description.setText(" ");
							Description.append("We now create a LJ Molecule\n");
							Description.append("With a 2-Centred-With-Quad Potential\n");
							EditVariables.setEnabled(true);
						}
					}
					if(PotentialType == "CO2"){
						
						potentialList = CO2PotentialClassList;
						EditVariables.setEnabled(false);
						table1.setEnabled(true);
						table2.setEnabled(true);
						
						if(list.getSelectedIndex() == 0){
							Description.setText(" ");
							Description.append("We now create a single CO2 Molecule\n");
							Description.append("With a 2-Centred-With-Quad Potential\n");
							
							}
					
						if(list.getSelectedIndex() == 1){
							Description.setText(" ");
							Description.append("We now create a single CO2 Molecule\n");
							Description.append("Modeled as EPM2\n");
							
							}
						if(list.getSelectedIndex() == 2){
							Description.setText(" ");
							Description.append("We now create a single CO2 Molecule\n");
							Description.append("Modeled as TRAPPE\n");
							}
					}
					if(PotentialType == "Methane"){
						potentialList = MethanePotentialClassList;
						EditVariables.setEnabled(false);
						
						
						table1.setEnabled(true);
						table2.setEnabled(true);
						obj[0] = new Integer(1);
						if(list.getSelectedIndex() == 0){
							Description.setText(" ");
							Description.append("We now create a single Methane Molecule\n");
							Description.append("Modelled as TRAPPE\n");
							}
						
					}
					
					if(PotentialType == "Ethane"){
						potentialList = EthanePotentialClassList;
						EditVariables.setEnabled(false);
						
						table1.setEnabled(true);
						table2.setEnabled(true);
						obj[0] = new Integer(2);
						if(list.getSelectedIndex() == 0){
							Description.setText(" ");
							Description.append("We now create a single Ethane Molecule\n");
							Description.append("Modelled as TRAPPE\n");
							}
						
					}
					
					if(PotentialType == "Propane"){
						potentialList = PropanePotentialClassList;
						EditVariables.setEnabled(false);
						
						table1.setEnabled(true);
						table2.setEnabled(true);
						obj[0] = new Integer(3);
						if(list.getSelectedIndex() == 0){
							Description.setText(" ");
							Description.append("We now create a single Propane Molecule\n");
							Description.append("Modelled as TRAPPE\n");
							}
						
					}
					if(PotentialType == "n-Alkane"){
						potentialList = nAlkanePotentialClassList;
						EditVariables.setEnabled(false);
						table1.setEnabled(true);
						table2.setEnabled(true);
						obj[0] = new Integer(4);
						if(list.getSelectedIndex() == 0){
							Description.setText(" ");
							Description.append("We now create a single n Molecule\n");
							Description.append("Modelled as TRAPPE\n");
						}
						if(!AddFlag){
							NAlkaneFrame1 = new JFrame("n-Alkane SpeciesA");
							if(NSpheres1 == null){
								NSpheres1 = new AlkaneSpheresUIView();
								NAlkaneFrame1.add(NSpheres1);
								NAlkaneFrame1.setMinimumSize(new Dimension(500,200));
								NAlkaneFrame1.setBounds(this.getParent().getParent().getWidth(), 0, 0, 0);
								NAlkaneFrame1.setVisible(true);
								NSpheres1.getCloseWindow().addActionListener(this);
								NSpheres1.getSaveValues().addActionListener(this);
							}
						}else{
							NAlkaneFrame2 = new JFrame("n-Alkane SpeciesB");
							if(NSpheres2 == null){
								NSpheres2 = new AlkaneSpheresUIView();
								NAlkaneFrame2.add(NSpheres2);
								NAlkaneFrame2.setMinimumSize(new Dimension(500,200));
								NAlkaneFrame2.setBounds(this.getParent().getParent().getWidth(), 0, 0, 0);
								NAlkaneFrame2.setVisible(true);
								NSpheres2.getCloseWindow().addActionListener(this);
								NSpheres2.getSaveValues().addActionListener(this);
							}
						}
						
						
					}
					
					if(PotentialType == "Water"){
						potentialList = WaterPotentialClassList;
						EditVariables.setEnabled(false);
						
						table1.setEnabled(true);
						table2.setEnabled(true);
						
						if(list.getSelectedIndex() == 0){
							Description.setText(" ");
							Description.append("We now create a single Water Molecule\n");
							Description.append("Modelled as SPC\n");
							}
						
					}
					
					if(!Remove.isEnabled()){
						Remove.setEnabled(true);
					}
					if(!AddSameSpecies.isEnabled()){
						AddSameSpecies.setEnabled(true);
					}
					if(!AddAnotherSpecies.isEnabled()){
						AddAnotherSpecies.setEnabled(true);
					}
					try{
						if(!AddFlag){
							TotalMoleculeCount = 0;
							if(potential1 == null){
								@SuppressWarnings("rawtypes")
								Constructor[] C = potentialList[list.getSelectedIndex()].getConstructors();
								if(C.length > 1){
									for(int i = 0; i < C.length;i++){
										if(C[i].getParameterTypes().length > 0){
											potential1 = (CreateSpeciesDM_IFactory) potentialList[list.getSelectedIndex()].getConstructor(Integer.TYPE).newInstance(obj);
										}
									}
								}
								else{
									potential1 = (CreateSpeciesDM_IFactory) potentialList[list.getSelectedIndex()].getConstructor().newInstance(new Object[0]);
								}
							}
							else{
								potential1 = null;
								/*for(int i = 0;i < 8;i++){
									if(sList.getObject(i) != null && AddFlag == false){
										sList.removeSpeciesAtIndex(i);
									}
								}*/
								if(sList.size() >= 1 && AddFlag == false){
									sList.clear();
								}
								
								//potential1 = (ParameterMapping) potentialList[list.getSelectedIndex()].getConstructor(Integer.TYPE).newInstance(new Object[0]);
								@SuppressWarnings("rawtypes")
								Constructor[] C = potentialList[list.getSelectedIndex()].getConstructors();
								if(C.length > 1){
									for(int i = 0; i < C.length;i++){
										if(C[i].getParameterTypes().length > 0){
											potential1 = (CreateSpeciesDM_IFactory) potentialList[list.getSelectedIndex()].getConstructor(Integer.TYPE).newInstance(obj);
										}
									}
								}
								else{
									potential1 = (CreateSpeciesDM_IFactory) potentialList[list.getSelectedIndex()].getConstructor().newInstance(new Object[0]);
								}
							
							}
							ParameterArray = potential1.getParametersArray();
							if(NSpheres1!= null){
								potential1.setParameter("NUMBER", Integer.toString(Alkane1Spheres));
							}
							//potential1.createSpeciesFactory();
							PotentialSite = potential1.getPotentialSites();
							ParamAndValue = potential1.getParamAndValues();
							
							
						}
						else{
							@SuppressWarnings("rawtypes")
							Constructor[] C = potentialList[list.getSelectedIndex()].getConstructors();
							if(C.length > 1){
								for(int i = 0; i < C.length;i++){
									if(C[i].getParameterTypes().length > 0){
										
											CreateSpeciesDM_IFactory TempPotential2 = (CreateSpeciesDM_IFactory) potentialList[list.getSelectedIndex()].getConstructor(Integer.TYPE).newInstance(obj);
											if(TempPotential2.getMoleculeDisplayName() == potential1.getMoleculeDisplayName() && NSpheres1 == null && NSpheres2 == null){
												//custom title, error icon
												MessageFrame = new JFrame("Alert");
												MessageAlert = new AlertMessageUIView(MessageFrame,"You ve chosen the same potential. Please choose another");
												MessageAlert.getCloseWindow().addActionListener(this);
												potential2 = null;
												return;
											}
											else{
												potential2 = (CreateSpeciesDM_IFactory) potentialList[list.getSelectedIndex()].getConstructor(Integer.TYPE).newInstance(obj);
											}
										
									}
								}
							}
							else{
								
									CreateSpeciesDM_IFactory TempPotential2 = (CreateSpeciesDM_IFactory) potentialList[list.getSelectedIndex()].getConstructor().newInstance(new Object[0]);
									if(TempPotential2.getMoleculeDisplayName() == potential1.getMoleculeDisplayName() && NSpheres1 == null && NSpheres2 == null){
										//custom title, error icon
										MessageFrame = new JFrame("Alert");
										MessageAlert = new AlertMessageUIView(MessageFrame,"You ve chosen the same potential. Please choose another");
										MessageAlert.getCloseWindow().addActionListener(this);
										return;
									}
									else{
										potential2 = (CreateSpeciesDM_IFactory) potentialList[list.getSelectedIndex()].getConstructor().newInstance(new Object[0]);
									}
								
							}
							ParameterArray = potential2.getParametersArray();
							if(NSpheres2!= null){
								potential2.setParameter("NUMBER", Integer.toString(Alkane2Spheres));
							}
							//potential2.createSpeciesFactory();
							PotentialSite = potential2.getPotentialSites();
							ParamAndValue = potential2.getParamAndValues();
							}
						}
					catch (Exception E){
						E.printStackTrace();
					}
					
					
					
					if (!AddFlag){
						//sList.addSpecies(potential1);
						sList.add(potential1.getMoleculeDisplayName());
						TotalMoleculeCount++;
						Molecule1Count++;
						
					}
					if (AddFlag){
						//sList.addSpecies(potential2);
						sList.add(potential2.getMoleculeDisplayName());
						TotalMoleculeCount++;
						Molecule2Count++;
					}
					
					displayMoleculeList();
					Description.append("1 Molecule added!\n");
					NoOfVariables = ParameterArray.length;
					if(!cellSelectionModel.isSelectionEmpty()){
						cellSelectionModel.clearSelection();
					}
					LJP1.removeData();
					for(int i = 0;i <ParamAndValue.length;i++){
						LJP1.setValueAt(ParamAndValue[i][0],i,0);
						LJP1.setValueAt(ParamAndValue[i][1],i,1);
					}
				}
				}	
			}
	
			
		
		
			if(e.getSource().equals(cellSelectionModel)){
				if( !e.getValueIsAdjusting()){
					
					String selectedParameter = null;
					int RowIndex = 0;
					if(table1.getSelectedRows().length != 0){
						int[] selectedRows = table1.getSelectedRows();
						RowIndex = selectedRows[0];
					}
					
					try{
							selectedParameter = (String)table1.getValueAt(RowIndex,0);
							for(PotentialParamDM_Description parameters : PotentialParamDM_Description.values()){
								if(selectedParameter.toUpperCase().contains(parameters.toString())){
									if(EditVariables.isEnabled()){
									Description.append("You have now chosen to edit the default value of " + parameters.toString() + "\n");}
									Description.append("This chosen parameter is " + parameters.Description()+"\n");
								}
							}
						}
						catch(NullPointerException n){}
					
					
					}
				
				}
		
			if(e.getSource().equals(cellSelectionModel2)){
				if( !e.getValueIsAdjusting()){
					
					String selectedParameter = null;
					String[][] ParamAndValue;
					int RowIndex = 0;
					if(table2.getSelectedRows().length != 0){
						int[] selectedRows = table2.getSelectedRows();
						RowIndex = selectedRows[0];
					}
					
					try{
							selectedParameter = (String)table2.getValueAt(RowIndex,1);
							if(selectedParameter == potential1.getMoleculeDisplayName()){
								ParamAndValue = potential1.getParamAndValues();
								
									
									for(int i = 0;i <ParamAndValue.length;i++){
										LJP1.setValueAt(ParamAndValue[i][0],i,0);
										LJP1.setValueAt(ParamAndValue[i][1],i,1);
									}
								
									
								
							}
							else{
								ParamAndValue = potential2.getParamAndValues();
								LJP1.removeData();
								for(int i = 0;i <ParamAndValue.length;i++){
									LJP1.setValueAt(ParamAndValue[i][0],i,0);
									LJP1.setValueAt(ParamAndValue[i][1],i,1);
								}
							}
						}
						catch(NullPointerException n){}
					
					
					}
				
			}

	}
	

	
	public void tableChanged(TableModelEvent e) {
		Object dataValue= null;
		Object dataString= null;
		int row = e.getFirstRow();
        int column = e.getColumn();
        TableModel model = (TableModel)e.getSource();
        
        if(model.equals(LJP1)){
        	try{
        	dataValue = model.getValueAt(row, column);
        	dataString = model.getValueAt(row, column - 1);
        	if(PotentialType == "LJ"){
        		if(PotentialJList.getSelectedIndex() == 0){	  
        			Description.append("You have now modified the default value!");
        			potential1.setParameter((String)dataString, (String)dataValue);
        		}
        	}
        	
        	}
        	catch(NullPointerException n){
        		
        	}
        	
        }
        if(model.equals(LJP2)){
        	
        }
        
	}
	
	
	
	//Listener method to run the simulation when the Run button is pressed
	void addSpeciesButtonListener(ActionListener mal) {
		AddSpecies.addActionListener(mal);
	}
	
	
	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				
				SpeciesSelectionUIView s = new SpeciesSelectionUIView();
				JFrame frame = new JFrame();
				frame.add(s);
				frame.setVisible(true);
				frame.setResizable(true);
				frame.setMinimumSize(new Dimension(750,700));
				frame.setBounds(0, 0,750, 700);
				
			}
		});
    }
}
