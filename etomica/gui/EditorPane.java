/**
 * EditorPane
 *
 * The EditorPane class is a splitpane that is the super pane of all component editor panes.  It is subclassed
 * by SimulationEditorPane and PotentialEditorPane.  The rightPane is a JList that holds all of the added
 * components of the current Simulation.instance field and the leftPane is a JPanel that holds JRadioButtons
 * corresponding to each available component for the SimulationEditorPane and JButtons corresponding to each
 * possible species-species interaction for the PotentialEditorPane.
 *
 * @author Bryan C. Mihalick
 * 1/19/01
 */
 
package etomica.gui;

import etomica.*;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JRadioButton;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;

public abstract class EditorPane extends javax.swing.JSplitPane {
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
     * Holds all constraints needed for displaying the next awt or swing component
     */
    final GridBagConstraints gbc = new GridBagConstraints();
    
    /**
     * Determines how to display an awt or swing component by using gbc from above
     */
    final GridBagLayout gbl = new GridBagLayout();
    
    /**
     * Button for adding objects to the JList of the right pane.  
     */
    final JButton addToSim = new JButton("Add");

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
     * Internal frame that lists the properties of a component
     */
    protected static PropertySheet propertySheet;
    
    /**
     * Gives the index of the currently selected object in the JList of the right pane
     */
    private int currentSelection;
        
    /**
     * If true, the property sheet is already visible I think
     */
    public static boolean added = false;
    
    /**
     * Buttongroup that provides mutual exclusion for the radio buttons of the left pane
     */
    final ButtonGroup simComponents = new ButtonGroup();
    
    /**
     * Handle to the simulation component that is going to be added to the JList on the right pane and
     * to the simulation.instance object
     */
    private Object component;
    
    /**
     * Name of editor pane
     */
    private String title;
    
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
    
    protected final SimulationEditor simulationEditor;

    public EditorPane(SimulationEditor ed){
        simulationEditor = ed;
    }
    
    public SimulationEditor simulationEditor() {return simulationEditor;}
    
    public int speciesCount(){ return simulationEditor().getSimulation().speciesCount(); }
    
    public javax.swing.DefaultListModel getComponentList(){ return componentList; }

    public void setEnabled(boolean b) {remove.setEnabled(b);}
    
    public void setTitle(String t){ 
        title = t; 
        if (title == "Phase")
            addToSim.setText("Add Phase");
    }
    public String getTitle() { return title; }

    public void setComponent(Object o) { component = o; }
    public Object getComponent() { return component; }
 
    public void setCurrentSelection(int cs) { currentSelection = cs; }
    public int getCurrentSelection() { return currentSelection; }
        
    /**
     * Provides a handle to the single instance of the property sheet so that it can be updated with the
     * new component's properties.
     */
    public static void setPropertySheet(PropertySheet p){ propertySheet = p; }
    
	protected class MyListSelectionListener implements javax.swing.event.ListSelectionListener, java.io.Serializable {
	    public void valueChanged(javax.swing.event.ListSelectionEvent lse){}
	}
	
	protected class MyInternalFrameAdapter extends javax.swing.event.InternalFrameAdapter implements java.io.Serializable {}

	protected class MyActionListener implements java.awt.event.ActionListener, java.io.Serializable {
	    public void actionPerformed(java.awt.event.ActionEvent ae){}
	}
}// end of EditorPane class