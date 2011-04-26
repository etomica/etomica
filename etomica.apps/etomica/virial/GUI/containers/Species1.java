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
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Locale;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.border.BevelBorder;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.swing.border.SoftBevelBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.plaf.metal.MetalLookAndFeel;
import javax.swing.plaf.metal.OceanTheme;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;

import com.jgoodies.looks.plastic.PlasticLookAndFeel;
import com.jgoodies.looks.plastic.PlasticXPLookAndFeel;
import com.jgoodies.looks.plastic.theme.ExperienceBlue;
import com.jgoodies.looks.plastic.theme.ExperienceGreen;

import etomica.virial.GUI.components.CreateP22CLJQ;
import etomica.virial.GUI.components.CreateP2LJ;
import etomica.virial.GUI.components.CreateP2LJQ;
import etomica.virial.GUI.components.SpeciesList;
import etomica.virial.GUI.models.LJParametersTableModel;
import etomica.virial.GUI.models.ParametersDouble;






public class Species1 extends JPanel implements SubPanelsInterface,ActionListener, ListSelectionListener, TableModelListener{
	
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
	
	//This needs to be removed once this code is ready to be added to the main code
	
	private JPanel TablePane;
	
	private String[] IntialList = {"---No species selected---","Press \'Add\' to select a species"};
	private String[] IntialPotentialList = {"---No potentials selected---","Press \'Add\' to select a potential"};
	
	private String[] SpeciesList = {"LJ","CO2","Methanol","Ethanol","Methane","Ethane","Propane","Naphthalene"};

	private String[] LJPotentialList = {"Spherical-2-Body","Spherical-2-Body-With-Q","2CLJQ"};
	private String[] CO2PotentialList = {"2CLJQ","TRAPPE-UnitedAtom","EPM2"};
	
	//For each species...
	private String PotentialType = null;
	private int NoOfVariables;
	private DialogBoxPanel DialogBox;
	
	//Parameter Tables
	private JTable table1;
	private JTable table2;
	
	//Description
	private JTextArea Description; 
	
	
	private JList SpeciesJList;
	private JList PotentialJList;
	
	private JButton ClearButton1;
	private JButton ClearButton2;
	private JButton ClearButton3;
	
	
	
	private LJParametersTableModel LJP1;
	private LJParametersTableModel LJP2;
	
	private ListSelectionModel cellSelectionModel;
	
	private CreateP2LJ p2LJ;
	private CreateP2LJQ p2LJQ;
	private CreateP22CLJQ p22CLJQ;
	
	private SpeciesList sList;
	private int ParameterListenerCallCount = 0;
	
	Species1(){
	super();
	initComponents();
	}
	public void initComponents(){
	this.setLayout(new BorderLayout());
	//
	System.out.println(this.getBounds().width);
	System.out.println(this.getHeight());
	
	
	JPanel panel1 = new JPanel();
	GridBagLayout gridbaglayout = new GridBagLayout();
	panel1.setLayout(gridbaglayout);
	this.add(panel1,BorderLayout.NORTH);
	
	AddSpecies = new JButton("Add");
	AddSpecies.setFocusPainted(false);
	AddSpecies.addActionListener(this);
	
	ClearButton1 = new JButton(new ImageIcon("Button-Close-icon.jpg"));
	ClearButton1.setPreferredSize(new Dimension(50,50));
	
	SpeciesJList = new JList(IntialList);
	SpeciesJList.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
	SpeciesJList.setLayoutOrientation(JList.HORIZONTAL_WRAP);
	SpeciesJList.setVisibleRowCount(-1);
	
 
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
   
    ClearButton2 = new JButton(new ImageIcon("Button-Close-icon.png"));
    ClearButton2.setPreferredSize(new Dimension(50,50));

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
    
    ClearButton3 = new JButton(new ImageIcon("Button-Close-icon.png"));
    ClearButton3.setPreferredSize(new Dimension(50,50));
    
    JPanel ParametersPane = new JPanel();
    ParametersPane.setLayout(new BoxLayout(ParametersPane, BoxLayout.PAGE_AXIS));
    
    TablePane = new JPanel();
    TablePane.setBorder(BorderFactory.createLoweredBevelBorder());
    
    LJP1 = new LJParametersTableModel();
    LJP2 = new LJParametersTableModel();
    
    //table1 = new JTable(LJP1.getData(),LJP1.columnNames);
    table1 = new JTable();
    table1.setModel(LJP1);
    table1.setEnabled(false);
   // table1.setCellSelectionEnabled(true);
    cellSelectionModel = table1.getSelectionModel();
    cellSelectionModel.addListSelectionListener(this);
    
    
    TableColumn column1 = null;
    for (int i = 0; i < 2; i++) {
        column1 = table1.getColumnModel().getColumn(i);
        column1.setPreferredWidth(100);
    }
    
    //table2 = new JTable(LJP2.getData(),LJP2.columnNames);
    table2 = new JTable();
    table2.setModel(LJP2);
    table2.setEnabled(false);
    
    
    TableColumn column2 = null;
    for (int i = 0; i < 2; i++) {
        column2 = table2.getColumnModel().getColumn(i);
        column2.setPreferredWidth(100);
       
    }
    //table1.setBounds(x, y, width, height)
    TablePane.add(table1);
    TablePane.add(table2);
    
    
    
    JLabel label3 = new JLabel("Parameters");
    label3.setForeground(Color.WHITE);
    JPanel labelPanel3 = new JPanel(new FlowLayout(FlowLayout.LEFT));
    labelPanel3.setBackground(new Color(0,78,152));
    labelPanel3.setBounds(ParametersPane.getBounds().x, ParametersPane.getBounds().y, 400, label3.getHeight());
    labelPanel3.add(label3);
    
    ParametersPane.add(labelPanel3);
    ParametersPane.add(Box.createRigidArea(new Dimension(5,0)));
    ParametersPane.add(TablePane);
    
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
    DescriptionButtonPane.add(AddSameSpecies);
    DescriptionButtonPane.add(AddAnotherSpecies);
    DescriptionButtonPane.add(Remove);
    
	
    JComponent[] componentLeft = {AddSpecies,AddPotential,EditVariables,new JPanel()};
    JComponent[] componentRight = {listPane,list2Pane,ParametersPane,DescriptionButtonPane};
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

	public String[] getLJPotentialList() {
		return LJPotentialList;
	}

	public void setLJPotentialList(String[] lJPotentialList) {
		LJPotentialList = lJPotentialList;
	}

	public String[] getCO2PotentialList() {
		return CO2PotentialList;
	}

	public void setCO2PotentialList(String[] cO2PotentialList) {
		CO2PotentialList = cO2PotentialList;
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

	public LJParametersTableModel getLJP1() {
		return LJP1;
	}

	public void setLJP1(LJParametersTableModel lJP1) {
		LJP1 = lJP1;
	}

	public LJParametersTableModel getLJP2() {
		return LJP2;
	}

	public void setLJP2(LJParametersTableModel lJP2) {
		LJP2 = lJP2;
	}

	public JButton getAddSpecies() {
		return AddSpecies;
	}

	public JButton getAddPotential() {
		return AddPotential;
	}

	public JButton getEditVariables() {
		return EditVariables;
	}

	/*
	
	public void addPressSpeciesButton(ActionListener B){
		AddSpecies.addActionListener(B);
	}
	
	public void addPressPotentialButton(ActionListener B){
		AddPotential.addActionListener(B);
	}
	
	public void addPressEditPotentialButton(ActionListener B){
		EditVariables.addActionListener(B);
	}*/
	
	public void actionPerformed(ActionEvent e){
		
		if(e.getSource().equals(AddSpecies)){
			SpeciesJList.setListData(SpeciesList);
			PotentialJList.setEnabled(true);
			AddPotential.setEnabled(true);
			sList = new SpeciesList();
		}
		
		if(e.getSource().equals(AddPotential)){
			
			if(SpeciesJList.getSelectedIndex() == 0){
				Description.setText(" ");
				Description.append("We now create a single LJ Molecule\n");
				PotentialJList.setListData(LJPotentialList);
				PotentialType = "LJ";
			}
			
			if(SpeciesJList.getSelectedIndex() == 1){
				Description.setText(" ");
				Description.append("We now create a single CO2 Molecule\n");
				PotentialJList.setListData(CO2PotentialList);
				PotentialType = "CO2";
			}
			
		}
		
		
		if(e.getSource().equals(EditVariables)){
			table1.setEnabled(true);
			LJP1.addTableModelListener(this);
			table2.setEnabled(true);
			LJP2.addTableModelListener(this);
			
		}
		
		if(e.getSource().equals(Remove)){
			if(sList.getId() != 0){
				sList.removeSpecies();
				Description.append("You ve just removed the last added species!!\n");
			}
			

		}
		
		if(e.getSource().equals(AddSameSpecies)){
			if(sList.getId() < 8){
				if(sList.getObject(sList.getId()-1) instanceof CreateP2LJ){
					CreateP2LJ clone = (CreateP2LJ)p2LJ.clone();
					sList.addSpecies(clone);
					
				}
				
				if(sList.getObject(sList.getId()-1) instanceof CreateP2LJQ){
					CreateP2LJQ clone = (CreateP2LJQ)p2LJQ.clone();
					sList.addSpecies(clone);
				}
				
				if(sList.getObject(sList.getId()-1) instanceof CreateP22CLJQ){
					CreateP22CLJQ clone = (CreateP22CLJQ)p22CLJQ.clone();
					sList.addSpecies(clone);
				}
			}
			Description.append("So far " + Integer.toString(sList.getId())+ " molecules have been added\n");
		}
	}
	
	public void valueChanged(ListSelectionEvent e){
		
			if(e.getSource().equals(PotentialJList)){
				if( !e.getValueIsAdjusting()){
				if(PotentialType == "LJ"){
					JList list = (JList)e.getSource();
					EditVariables.setEnabled(true);
					String[] ParameterArray = null;
				
				
				
					if(list.getSelectedIndex() == 0){
						Description.setText(" ");
						Description.append("We now create a single LJ Molecule\n");
						Description.append("With a Spherical-2-Body Potential\n");
						try{
							if(p2LJ == null){
								p2LJ = new CreateP2LJ();
								p2LJ.createP2LJSpeciesFactory();
								
							}
							ParameterArray = p2LJ.getParametersArray();
						}catch(Exception E){
						}finally{
							if(p2LJQ != null){
								p2LJQ = null;
								sList.removeSpecies();
							}
							if(p22CLJQ != null){
								p22CLJQ = null;
								sList.removeSpecies();
							}
							sList.addSpecies(p2LJ);
						}
					}
				
					if(list.getSelectedIndex() == 1){
						Description.setText(" ");
						Description.append("We now create a single LJ Molecule\n");
						Description.append("With a Spherical-2-Body-With-Quad Potential\n");
						if(p2LJQ == null){
							
							if(p2LJ != null){
								
								p2LJ = null;
								sList.removeSpecies();
							}
							if(p22CLJQ != null){
								p22CLJQ = null;
								sList.removeSpecies();
							}
							p2LJQ = new CreateP2LJQ();
							p2LJQ.createP2LJQSpecies();
							sList.addSpecies(p2LJQ);
						}
						ParameterArray = p2LJQ.getParametersArray();
					}
				
					if(list.getSelectedIndex() == 2){
					
						Description.setText(" ");
						Description.append("We now create a LJ Molecule\n");
						Description.append("With a 2-Centred-With-Quad Potential\n");
						if(p22CLJQ == null){
							if(p2LJ != null){
								
								p2LJ = null;
								sList.removeSpecies();
							}
							if(p2LJQ != null){
								p2LJQ = null;
								sList.removeSpecies();
							}
							p2LJQ = new CreateP2LJQ();
							p2LJQ.createP2LJQSpecies();
							sList.addSpecies(p2LJQ);
						}
					//LJP.UpdateObjectData(p22CLJQ.getParametersArray());
						ParameterArray = p22CLJQ.getParametersArray();
					
					}
					NoOfVariables = ParameterArray.length;
					LJP1.removeData();
					for(int i = 0;i < NoOfVariables; i++){
						for(ParametersDouble parameters : ParametersDouble.values()){
							if(ParameterArray[i]==parameters.toString()){
								LJP1.setValueAt(parameters.toString().toLowerCase(Locale.ENGLISH), i, 0);
								LJP1.setValueAt(Double.toString(parameters.DefaultValue()), i, 1);
							}
						}
					}
					
				}
			
				if(PotentialType == "CO2"){
				}
			
				if(PotentialType == "Alkane"){
				
				}
				if(PotentialType == "Alcohol"){
				}
				}
			}
		
		
			if(e.getSource().equals(cellSelectionModel)){
				if( !e.getValueIsAdjusting()){
					String selectedParameter = null;
					int[] selectedRows = table1.getSelectedRows();
					int RowIndex = selectedRows[0];
					try{
						selectedParameter = (String)table1.getValueAt(RowIndex,0);
						for(ParametersDouble parameters : ParametersDouble.values()){
							if(parameters.toString().equals(selectedParameter.toUpperCase())){
								Description.append("You have now chosen to edit the default value of " + parameters.toString() + "\n");
								Description.append("This chosen parameter is " + parameters.Description()+"\n");
							}
						}
					}
					catch(NullPointerException n){}
					}  	
			}
			
		
		
		
		
	}
	
	@Override
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
        			p2LJ.setParameter((String)dataString, (String)dataValue);
        		}
        	}
        	}
        	catch(NullPointerException n){
        		
        	}
        	
        }
        if(model.equals(LJP2)){
        	
        }
        
	}

	
	
	
	
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				
				Species1 s = new Species1();
				JFrame frame = new JFrame();
				frame.add(s);
				frame.setVisible(true);
				frame.setResizable(true);
				frame.setMinimumSize(new Dimension(530,600));
				frame.setBounds(0, 0, 530, 600);
				
			}
		});
    }

	

	

}
