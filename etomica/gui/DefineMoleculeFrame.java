/**
 * DefineMoleculeFrame
 *
 * The DefineMoleculeFrame class is responsible for creating a new JSplitPane that gives the user
 * choices for creating a new Molecule to add to the simulation.
 *
 * @author Bryan C. Mihalick
 * 9/18/00
 */
 
package etomica.gui;

import etomica.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import javax.swing.JButton;
import java.io.File;
import java.io.FilenameFilter;

public class DefineMoleculeFrame extends javax.swing.JInternalFrame {
    /**
     * Static array of all simulation components that extend AtomType.class
     */
    public static Class[] atomTypeClasses;
    
    SimulationEditor simulationEditor;
        
    DefineMoleculeFrame(SimulationEditor ed){
        simulationEditor = ed;
        setResizable(true);
        setMaximizable(true);
        setIconifiable(true);
        setClosable(false);
        setVisible(true);
        setBounds(515, 60, 320, 200);
        setTitle("Define New Molecule");
        
        getContentPane().add(new DefineMoleculePane(this));
    }
        
    class DefineMoleculePane extends javax.swing.JSplitPane {
        /**
         * Determines format of the JList on the left pane of the DefineMolecule SplitPane
         */
        final javax.swing.DefaultListModel atomTypeList = new javax.swing.DefaultListModel();    

        /**
         * Determines format of the JList on the right pane of the DefineMolecule SplitPane
         */
        final javax.swing.DefaultListModel atomTypeArrayList = new javax.swing.DefaultListModel();    

        /**
         * JList that contains all of the simulation components that were added from the left pane's choices
         */
        final javax.swing.JList leftPaneList = new javax.swing.JList(atomTypeList);
        
        /**
         * JList that contains all of the simulation components that were added from the left pane's choices
         */
        final javax.swing.JList rightPaneList = new javax.swing.JList(atomTypeArrayList);
        
        /**
         * Lists all of the simulation components corresponding to the respective tabs name.  These are listed
         * as radio buttons.
         */
        final javax.swing.JScrollPane leftPane = new javax.swing.JScrollPane(leftPaneList);
        
        /**
         * Scrollable pane that holds the rightPaneList from above
         */
        final javax.swing.JScrollPane rightPane = new javax.swing.JScrollPane(rightPaneList);
        
        /**
         * JPanel that contains the rightPane and the Add, Remove, and setConfiguration buttons.
         */
        final javax.swing.JPanel rightPanel = new javax.swing.JPanel();
        
        /**
        * Holds all constraints needed for displaying the next awt or swing component
        */
        final GridBagConstraints gbc = new GridBagConstraints();
        
        /**
        * Determines how to display an awt or swing component by using gbc from above
        */
        final GridBagLayout gbl = new GridBagLayout();
       
        /**
         * Gives the index of the currently selected object in the AtomType JList of the left pane
         */
        int currentAtomType;
     
        /**
         * Gives the index of the currently selected object in the AtomTypeArray JList of the right pane
         */
        int currentAtomTypeArray;
        
        /**
         * Holds AtomTypes for exportation to the Configuration Panel
         */
        AtomType types[];
        
        /**
         * Adds currently selected AtomType to the AtomTypeArray
         */
        JButton addAtomType;
        
        /**
         * Removes currently selected AtomType from the AtomTypeArray
         */
        JButton removeAtomType;
        
        /**
         * Opens DefineConfigurationFrame for easy configuring of newly created molecule
         */
        JButton setConfiguration;
        
        /**
         * Constructor that sets up all the properties of the JSplitPane including the AtomType list(holds
         * all possible choices of AtomType) and the AtomType Array List(holds all AtomTypes chosen from the
         * AtomType list).
         */
        public DefineMoleculePane(javax.swing.JInternalFrame frame){
            final javax.swing.JInternalFrame parentFrame = frame;
            setBounds(515, 60, 350, 200);
    		
    		
    		leftPane.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
    		rightPane.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
		    leftPaneList.getSelectionModel().setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
		    leftPaneList.addListSelectionListener(new javax.swing.event.ListSelectionListener(){
		        public void valueChanged(javax.swing.event.ListSelectionEvent lse){
		            currentAtomType = leftPaneList.getLeadSelectionIndex();
		        }});
		    rightPaneList.getSelectionModel().setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
		    rightPaneList.addListSelectionListener(new javax.swing.event.ListSelectionListener(){
		        public void valueChanged(javax.swing.event.ListSelectionEvent lse){
                    currentAtomTypeArray = rightPaneList.getLeadSelectionIndex();
		        }});
   		
   		    for (int i = 0; i < atomTypeClasses.length; i++) {
   		        try {
		            atomTypeList.addElement(((Class)atomTypeClasses[i]).newInstance());
                }
	            catch(InstantiationException exc) {}
	            catch(IllegalAccessException exc) {}                
            }
		    atomTypeList.addElement(new AtomType.Well(1.0, java.awt.Color.red,1.0,0.1));
		    atomTypeList.addElement(new AtomType.Wall(1.0, java.awt.Color.red,0.0,0));
    		
    		rightPanel.setLayout(gbl);
    		gbc.gridx = 0;
    		gbc.gridy = 0;
    		gbc.weightx = 95;
    		gbc.weighty = 95;
    		gbc.gridwidth = 4;
    		gbc.gridheight = 1;
            gbc.fill = gbc.BOTH;
            gbc.anchor = gbc.CENTER;
    		gbl.addLayoutComponent(rightPane, gbc);
    		rightPanel.add(rightPane);
    		
            addAtomType = new JButton("Add");
            addAtomType.addActionListener(new ActionListener(){
                public void actionPerformed(ActionEvent e){
                    atomTypeArrayList.addElement(atomTypeList.getElementAt(currentAtomType));
                    removeAtomType.setEnabled(true);
                }});
            gbc.gridx++;
            gbc.gridy++;
            gbc.gridwidth = 1;
            gbc.weightx = 0;
            gbc.weighty = 0;
            gbc.fill = gbc.NONE;
    		gbl.addLayoutComponent(addAtomType, gbc);
            rightPanel.add(addAtomType);
            
            removeAtomType = new JButton("Remove");
            removeAtomType.addActionListener(new ActionListener(){
                public void actionPerformed(ActionEvent e){
                    if (currentAtomTypeArray != -1){
                        atomTypeArrayList.removeElementAt(currentAtomTypeArray);
                        currentAtomTypeArray = -1;
                        if (atomTypeArrayList.getSize() == 0){
                            removeAtomType.setEnabled(false);
                        }
                    }
                }});
    		gbc.gridx += 2;
    		gbl.addLayoutComponent(removeAtomType, gbc);
            rightPanel.add(removeAtomType);
            
            setConfiguration = new JButton("Set Configuration");
            setConfiguration.addActionListener(new ActionListener(){
                public void actionPerformed(ActionEvent e){
                    types = new AtomType[atomTypeArrayList.size()];
                    int i = 0;
                    for (java.util.Enumeration enum = atomTypeArrayList.elements(); enum.hasMoreElements();){
                        types[i++] = ((AtomType)enum.nextElement());
                    }
                    try {
                        parentFrame.setClosed(true);
                    }
                    catch (java.beans.PropertyVetoException pve) {}
                    DefineConfigurationFrame configFrame = new DefineConfigurationFrame(simulationEditor);
                    configFrame.setAtomTypes(types);
                    Etomica.DesktopFrame.desktop.add(configFrame);
                    try {
                        configFrame.setSelected(true);
                    }
                    catch (java.beans.PropertyVetoException pve) {}
                }});
            gbc.gridx = 0;
            gbc.gridy++;
            gbc.gridwidth = 4;
    		gbl.addLayoutComponent(setConfiguration, gbc);
            rightPanel.add(setConfiguration);
            
            leftPane.setMinimumSize(new java.awt.Dimension(150, 80));
		    setLeftComponent(leftPane);
            setRightComponent(rightPanel);
            this.setDividerLocation(0.5);            
        }// end of DefineMoleculePane constructor
    }// end of DefineMoleculePane class
    
    static {
        // Initialization of spaceClasses array
	    File dir = new File(etomica.Default.CLASS_DIRECTORY);
	    String[] files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
	            return name.startsWith("AtomType")
	                && name.endsWith("class")
	                && !name.startsWith("AtomType$Rod")
	                && !name.startsWith("AtomType$Disk")
	                && !name.equals("AtomType$Carbon.class")
	                && !name.equals("AtomType.class");}
	        });
	    atomTypeClasses = new Class[files.length];
	    for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        atomTypeClasses[i] = null;
	        try{
	            atomTypeClasses[i] = Class.forName("etomica."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End of initialization of atomTypeClasses array
	}
}// end of DefineMoleculeFrame class