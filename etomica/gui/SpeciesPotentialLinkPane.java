/**
 * SpeciesPotentialLinkPane
 *
 * The SpeciesPotentialLinkPane class is a splitpane that lists all the simulation components of a 
 * respective potential category (potential 1 or potential 2) on the leftside, and all of the added 
 * components from the leftside list on the rightside in a JList format.
 *
 * @author Bryan C. Mihalick
 * 8/16/00
 */

package simulate.gui;

import simulate.*;
import java.awt.Color;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JRadioButton;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import java.io.File;
import java.io.FilenameFilter;

public abstract class SpeciesPotentialLinkPane extends javax.swing.JSplitPane implements SimulationEditor.EditorPane {
    /**
     * Determines format of the JList on the right pane of the SimulationEditorPane
     */
    protected final javax.swing.DefaultListModel componentList = new javax.swing.DefaultListModel();    
    
    /**
     * JList that contains all of the simulation components that were added from the left pane's choices
     */
    final javax.swing.JList rightPaneList = new javax.swing.JList(componentList);
    
    /**
     * Lists all of the simulation components corresponding to the respective tabs name.  These are listed
     * as radio buttons.
     */
    final javax.swing.JPanel leftPanePanel = new javax.swing.JPanel();
    
    /**
     * Scrollable pane that holds the leftPanePanel from above
     */
    final javax.swing.JScrollPane leftPane = new javax.swing.JScrollPane(leftPanePanel);
    
    /**
     * Scrollable pane that holds the rightPaneList from above
     */
    final javax.swing.JScrollPane rightPane = new javax.swing.JScrollPane(rightPaneList);

    /**
     * Buttongroup that provides mutual exclusion for the radio buttons of the left pane
     */
    final ButtonGroup simComponents = new ButtonGroup();
    
    /**
     * Makes it possible to determine which radio button was selected when the "add" button is pressed
     */
    final ButtonListener buttonListener = new ButtonListener();
    
    /**
     * Holds all constraints needed for displaying the next awt or swing component
     */
    final GridBagConstraints gbc = new GridBagConstraints();
    
    /**
     * Determines how to display an awt or swing component by using gbc from above
     */
    final GridBagLayout gbl = new GridBagLayout();

    /**
     * Handle to the simulation component that is going to be added to the JList on the right pane and
     * to the simulation.instance object
     */
    Object component;
    
    /**
     * Number of species added from the species pane
     */
    private int numOSpecies = 0;
    
    /**
     * Button for removing objects from the JList of the right pane.  Only enabled if an object is in the
     * JList
     */
    final JButton remove = new JButton("Remove");
    
    /**
     * Button for starting the simulation.  Only enabled if a sufficient number of simulation components
     * have been added to make a working simulation
     */
    final JButton start = new JButton("Start");
    
    /**
     * Array of all current JButtons that are pressed 
     */
    public static JButton[] currentButtons = new JButton[1];
    
    /**
     * Array of all potential JButtons on the left pane
     */
    public static JButton[] potButtons = new JButton[1];
    
    /**
     * Number of buttons on left pane that are currently pressed
     */
    int buttonCount = 0;
    
    /**
     * Index of the potential class that the currently selected JButton in the left pane corresponds to
     */
    private int currentIndex = 0;
    
    /**
     * Index of object that is currently selected in the JList on the right pane
     */
    private int currentSelection = -1;
    
    /**
     * The title that will be displayed on the property sheet when a component is selected
     */
    private String title;
    
    /**
     * If true, the current simulation component has already been added
     */
    public static boolean added = false;
    
    /**
     * Frame containing radioButtons for all of the Potential classes
     */
    public PotentialFrame potentialFrame;
    
    /**
     * If true, at least one species was added to the simulation.  This makes sure that the
     * actionListener of the remove button is only added once.
     */
    boolean firstSpecies = true;
    
    /**
     * Envelopes a simulation component in order to have the component's properties listed in the 
     * property sheet.
     */
    protected static Wrapper wrapper = null;
    
    /**
     * Internal frame that lists the properties of a component
     */
    protected static PropertySheet propertySheet;
    
    /**
     * Determines if Potential Arrays need to be instantiated or not
     */
    public static boolean makePotArrays = true;
    
    /**
     * Static array of all simulation components that extend potential1.class
     */
    public static Class[] potential1Classes;
    
    /**
     * Static array of all simulation components that extend potential2.class
     */
    public static Class[] potential2Classes;
        
    /**
     * Initial height of splitPane
     */
    protected int splitPaneHeight = 580;
    
    /**
     * Initial width of splitPane
     */
    protected int splitPaneWidth = 580;
    
    /**
     * Height of leftPane.  Used for determining when scrollbars should be added to the scrollPane
     */
    protected int leftPanelHeight;
    
    /**
     * Width of leftPane.  Used for determining when scrollbars should be added to the scrollPane
     */
    protected int leftPanelWidth;
    
    /**
     * Height of one JRadioButton.  Used for determining leftPanelHeight.
     */
    static final protected int radioButtonHeight = 24;

    /**
     * Height of one JButton.  Used for determining leftPanelHeight.
     */
    static final protected int jButtonHeight = 32;

    protected SimulationEditor simulationEditor;
    /**
     * Constructor that creates all of the left pane's radiobuttons and JButtons, as well as, the right 
     * pane's scrollpane and JList.  It also creates all the listeners for these swing components so that
     * the simulation can be updated as needed.
     */
    public SpeciesPotentialLinkPane(SimulationEditor ed){
        simulationEditor = ed;
        setSize(580, 580);
		setLeftComponent(leftPane);
        setRightComponent(rightPane);
        setDividerLocation(0.65);
        leftPanelWidth = (int)(0.65*splitPaneWidth)-15;
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
        gbc.gridwidth = gbc.REMAINDER;
        gbc.gridx = 0;
        gbc.gridy = 0;
        
        /**
         * This section sets up the JList of the right pane.  When an object from the list is selected,
         * its properties are displayed in the property sheet. 
         */
        rightPaneList.addListSelectionListener(new MyListSelectionListener(){
            public void valueChanged(javax.swing.event.ListSelectionEvent lse){
		        currentSelection = rightPaneList.getLeadSelectionIndex();
                for (int i = 0; i < leftPanePanel.getComponentCount(); i++){
                    if (leftPanePanel.getComponent(i) instanceof SpeciesPairButton){
                        if (((SpeciesPairButton)leftPanePanel.getComponent(i)).index == currentSelection && currentSelection != -1)
                            ((javax.swing.AbstractButton)leftPanePanel.getComponent(i)).setBackground(Color.red);
                        else ((javax.swing.AbstractButton)leftPanePanel.getComponent(i)).setBackground(Color.lightGray);
                    }
                }
                if (added == false)
                    ViewActions.PROPERTYLIST.actionPerformed(new ActionEvent(this, 0, ""));
	            if (rightPaneList.getSelectedValue() instanceof P1DefinedPotential) {
	                PotentialViewer potentialViewer = new PotentialViewer(simulationEditor);
	                potentialViewer.setParam(((P1DefinedPotential)rightPaneList.getSelectedValue()).atomPairPotArray,propertySheet,((P1DefinedPotential)rightPaneList.getSelectedValue()).getName());    
                    potentialViewer.setTitle("PotentialViewer - " + potentialViewer.title);
                    potentialViewer.setVisible(true);
	            }
	            else if (rightPaneList.getSelectedValue() instanceof P2DefinedPotential) {
	                PotentialViewer potentialViewer = new PotentialViewer(simulationEditor);
	                potentialViewer.setParam(((P2DefinedPotential)rightPaneList.getSelectedValue()).atomPairPotArray,propertySheet,((P2DefinedPotential)rightPaneList.getSelectedValue()).getName());    
                    potentialViewer.setTitle("PotentialViewer - " + potentialViewer.title);
                    potentialViewer.setVisible(true);
	            }
	            else {
	                if (rightPaneList.getSelectedValue() != null) {
	                    /*if (rightPaneList.getSelectedValue() instanceof Potential2) 
                            wrapper = new Wrapper(((Potential2)rightPaneList.getSelectedValue()).getOnlyPotential(), rightPaneList.getSelectedValue().toString(), "simulate.gui." + rightPaneList.getSelectedValue().toString()); 
                        else if (rightPaneList.getSelectedValue() instanceof Potential1)
                            wrapper = new Wrapper(((Potential1)rightPaneList.getSelectedValue()).getOnlyPotential(), rightPaneList.getSelectedValue().toString(), "simulate.gui." + rightPaneList.getSelectedValue().toString()); 
                        else */wrapper = new Wrapper(rightPaneList.getSelectedValue(), rightPaneList.getSelectedValue().toString(), "simulate.gui." + rightPaneList.getSelectedValue().toString());
                        //wrapper = new Wrapper(rightPaneList.getSelectedValue(), title, "simulate.gui." + title); 
                        propertySheet.setTarget(wrapper);
                        try {
                            propertySheet.setSelected(true);
                        }
                        catch(java.beans.PropertyVetoException pve){}
                    }
	                propertySheet.addInternalFrameListener(new MyInternalFrameAdapter(){
	                    public void internalFrameClosed( InternalFrameEvent ife ){
	                        rightPaneList.clearSelection();
	                    }});
                }
		    }});// end JList instantiation and setup
		rightPaneList.getSelectionModel().setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
	    leftPanePanel.setMinimumSize(new java.awt.Dimension(300, 200));
    }

    public SimulationEditor simulationEditor() {return simulationEditor;}
    
    /**
     * This method is called when a species is added or removed from the species pane.  It figures out
     * the number of possible potential interactions, creates buttons for each of these interactions,
     * determines the necessary layout of the buttons, and displays them.
     */
    public abstract void update();

    public void addButtons() {
        /**
         * This section adds the start button, but it's disabled until enough simulation components
         * have been added to make a working simulation.
         */
        gbc.gridx = 0;
        gbc.gridy++;
	    start.setEnabled(false);
	    start.addActionListener(new MyActionListener(){
	        public void actionPerformed(ActionEvent e){
	            SimulationFrame frame = (SimulationFrame)Etomica.simulationFrames.get(simulationEditor.getSimulation());
	            frame.setVisible(true);
                try { frame.setSelected(true); }
                catch (java.beans.PropertyVetoException exc){} // attempt was vetoed
            }});
        gbl.setConstraints(start,gbc);
        leftPanePanel.add(start);
        // end of start button
            
        /**
         * This section adds the remove button, but it's disabled until an object is added to the JList
         */
        gbc.gridx++;
        gbc.gridwidth = 2;
	    remove.setEnabled(false);
	    if (componentList.getSize() != 0){
	        remove.setEnabled(true);
	    }
	    if (firstSpecies == true) {
	        firstSpecies = false;
	        remove.addActionListener(new MyActionListener(){
	            public void actionPerformed(ActionEvent e){
	                // Checks what buttons corresponded to the deleted potential and resets those
	                // buttons so that a new potential can be designated for them.
	                for (int i = 0; i < potButtons.length-1; i++){
	                    if (((SpeciesPairButton)potButtons[i]).index == currentSelection){
	                        potButtons[i].setEnabled(true);
	                        potButtons[i].setBackground(Color.lightGray);
	                        ((SpeciesPairButton)potButtons[i]).index = -1;
	                        ((SpeciesPairButton)potButtons[i]).potential = null;
	                    }
	                }
	                currentIndex--;                         // Decreases the total number of current potentials
	                shiftIndices(currentSelection);         // Updates indices of each button
	                if (title == "Potential1") {
	                        Simulation.instance.unregister((Potential1)componentList.getElementAt(currentSelection));
                    }
                    else {
                        Simulation.instance.unregister((Potential2)componentList.getElementAt(currentSelection));
	                }
	                componentList.remove(currentSelection); // Removes selected potential from the component list
                    wrapper = new Wrapper(FileActions.OPEN, "null", "null"); 
                    propertySheet.setTarget(wrapper);
                    if (componentList.getSize() == 0)
                        ((JButton)e.getSource()).setEnabled(false);
                        
                    // Check if a sufficient number of components are added to allow a working simulation.
                    // If so, enable the start button.
	                if (simulationEditor.allRemoveEnabled())
	                    simulationEditor.setAllStart(true);
	                else simulationEditor.setAllStart(false);
	            }});
	    }
	    gbl.setConstraints(remove,gbc);
	    leftPanePanel.add(remove);
    	// end of remove button
    	    
    	/**
    	 * This section adds the property sheet button, that when pressed opens the property sheet
    	 * with the currently selected component's properties listed.  If no component is selected,
    	 * it still opens but is blank
    	 */
	    gbc.gridx--;
	    gbc.gridy++;
	    gbc.gridwidth = 3;
	    JButton propSheet = new JButton("Property Sheet");
    	propSheet.addActionListener(new MyActionListener(){
	        public void actionPerformed(ActionEvent e){
                if (added == false)     // If a property list hasn't been made yet, it makes one
                    ViewActions.PROPERTYLIST.actionPerformed(new ActionEvent(this, 0, ""));
	            if (rightPaneList.getSelectedValue() != null){ 
                    wrapper = new Wrapper(rightPaneList.getSelectedValue(), title, "simulate.gui." + title); 
                    propertySheet.setTarget(wrapper);   // Updates property sheet with the selected object's properties
                }
                else {
                    wrapper = new Wrapper(FileActions.OPEN, "null", "null"); 
                    propertySheet.setTarget(wrapper);
                }
                try {
                    propertySheet.setSelected(true);
                }
                catch(java.beans.PropertyVetoException pve){}
                propertySheet.setVisible(true);
                propertySheet.repaint();
	        }});
	    gbl.setConstraints(propSheet,gbc);
	    leftPanePanel.add(propSheet);
	    // end of property sheet button
    }
    
    /**
     * This method shifts the index of each JButton when a potential is removed from the JList
     */ 
    public void shiftIndices(int deletedIndex){
    
        for (int i = 0; i < potButtons.length-1; i++){
            if (((SpeciesPairButton)potButtons[i]).index > deletedIndex)
                ((SpeciesPairButton)potButtons[i]).index--;
        }
    }// end of shiftIndices
    
    /**
     * Provides a handle to the single instance of the property sheet so that it can be updated with the
     * new component's properties.
     */
    public static void setPropertySheet(PropertySheet p){ propertySheet = p; }
    
    public void setNumOSpecies(int n){ numOSpecies = n; }
    public int getNumOSpecies(){ return numOSpecies; }
    
    public void setComponent(Object o) { component = o; }
    public Object getComponent() { return component; }
 
    public void setCurrentIndex(int cs) { currentIndex = cs; }
    public int getCurrentIndex() { return currentIndex; }
    
    public void setCurrentSelection(int cs) { currentSelection = cs; }
    public int getCurrentSelection() { return currentSelection; }

    public void setTitle(String t){ title = t; }
    public String getTitle() { return title; }
    
	public javax.swing.DefaultListModel getComponentList(){ return componentList; }
    
    /**
     * This class allows buttons to hold a handle to the potential class that is associated to the button,
     * as well as the corresponding index of that potential in the JList
     */
    public class SpeciesPairButton extends javax.swing.JButton{
        Class potential = null;
        int index = -1;
        int speciesIndex1 = -1;
        int speciesIndex2 = -1;
        
        SpeciesPairButton(String name){
            super(name);
        }// end of SpeciesPairButton constructor        
    }// end of SpeciesPairButton class
    
    /**
     * Extension class of ActionListener that saves a handle to the buttons that are currently 
     * selected.  It also turns the selected buttons to red, and then instantiates a PotentialFrame
     * window that lists all of the possible potential classes that can be added to the simulation.
     */
	protected class ButtonListener implements ActionListener, java.io.Serializable {
	    
	    public void actionPerformed(ActionEvent evt) {
	        currentButtons[buttonCount] = ((JButton)evt.getSource());
	        DefinePotentialFrame.setSpeciesIndex1(((SpeciesPairButton)currentButtons[buttonCount]).speciesIndex1);
	        DefinePotentialFrame.setSpecies1(((Species)simulationEditor.speciesEditor.componentList.getElementAt(((SpeciesPairButton)currentButtons[buttonCount]).speciesIndex1)));
	        DefinePotentialFrame.setSpecies2(((Species)simulationEditor.speciesEditor.componentList.getElementAt(((SpeciesPairButton)currentButtons[buttonCount]).speciesIndex1)));
            PotentialFrame.atomPairPotArray = new Class[DefinePotentialFrame.species1.getAtomsPerMolecule()][DefinePotentialFrame.species1.getAtomsPerMolecule()];
	        try {
	            DefinePotentialFrame.setSpeciesIndex2(((SpeciesPairButton)currentButtons[buttonCount]).speciesIndex2);
                PotentialFrame.atomPairPotArray = new Class[DefinePotentialFrame.species1.getAtomsPerMolecule()][DefinePotentialFrame.species2.getAtomsPerMolecule()];
	        }
	        catch (java.lang.ArrayIndexOutOfBoundsException exc) {}
	        
	        
	        ((SpeciesPairButton)currentButtons[buttonCount]).index = currentIndex;
	        buttonCount++;
            currentSelection = -1;
	        rightPaneList.clearSelection();
	        ((java.awt.Component)evt.getSource()).setBackground(Color.red);

            if (makePotArrays) {
///                Simulation.potential1 = new Potential1[getNumOSpecies()];
///                Simulation.potential2 = new Potential2[getNumOSpecies()][getNumOSpecies()];
                makePotArrays = false;
            }
            
	        if (buttonCount == 1){
	            PotentialFrame.setRemoveButton(remove);
	            if (getTitle() == "Potential1") {
	                potentialFrame = new PotentialFrame(simulationEditor);
	                potentialFrame.setTitle(getTitle());
	                potentialFrame.setPotentialEditor(simulationEditor.potential1Editor);
	            }
	            else {
	                potentialFrame = new PotentialFrame(simulationEditor);
	                potentialFrame.setTitle(getTitle());
	                potentialFrame.setPotentialEditor(simulationEditor.potential2Editor);
	            }
	            Etomica.DesktopFrame.desktop.add(potentialFrame);
	            try {
	                potentialFrame.setSelected(true);
	            }
	            catch(java.beans.PropertyVetoException pve){}
            }
        }// end of actionPerformed
	}// end of ButtonListener class
	
	private class MyListSelectionListener implements javax.swing.event.ListSelectionListener, java.io.Serializable {
	    public void valueChanged(javax.swing.event.ListSelectionEvent lse){}
	}
	
	public class MyInternalFrameAdapter extends javax.swing.event.InternalFrameAdapter implements java.io.Serializable {}

	public class MyActionListener implements java.awt.event.ActionListener, java.io.Serializable {
	    public void actionPerformed(java.awt.event.ActionEvent ae){}
	}

	static {
    	// Initialization of potential1Classes array
	    File dir = new File(simulate.Default.CLASS_DIRECTORY);
	    String[] files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("P1")
                    && name.endsWith("class")
	                && !name.startsWith("P1$");}
	        });
        potential1Classes = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        potential1Classes[i] = null;
	        try{
	            potential1Classes[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of potential1Classes array
    	
    	// Initialization of potential2Classes array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("P2")
	                && name.endsWith("class")
	                && !name.startsWith("P2$");}
	        });
        potential2Classes = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        potential2Classes[i] = null;
	        try{
	            potential2Classes[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of potential2Classes array
	}// End of static block
}// end of SpeciesPotentialLinkPane class