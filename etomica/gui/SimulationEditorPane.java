/**
 * SimulationEditorPane
 *
 * The SimulationEditorPane class is a splitpane that lists all the simulation components of a respective
 * category (species, integrator, phase, controller, device, display, meter, but not potential) on the 
 * leftside, and all of the added components from the leftside list on the rightside in a JList format.
 *
 * @author Bryan C. Mihalick
 * 8/16/00
 */
 
package simulate.gui;

import simulate.*;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JInternalFrame;
import javax.swing.JRadioButton;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public abstract class SimulationEditorPane extends javax.swing.JSplitPane {
    /**
     * Determines format of the JList on the right pane of the SimulationEditorPane
     */
    final javax.swing.DefaultListModel componentList = new javax.swing.DefaultListModel();    
    
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
     * Envelopes a simulation component in order to have the component's properties listed in the 
     * property sheet.
     */
    protected static Wrapper wrapper = null;
    
    /**
     * Internal frame that lists the properties of a component
     */
    protected static PropertySheet propertySheet;
    
    /**
     * Gives the index of the currently selected object in the JList of the right pane
     */
    private int currentSelection;
        
    /**
     * If true, the current simulation component has already been added
     */
    public static boolean added = false;
    
    /**
     * Buttongroup that provides mutual exclusion for the radio buttons of the left pane
     */
    final ButtonGroup simComponents = new ButtonGroup();
    
    /**
     * Makes it possible to determine which radio button was selected when the "add" button is pressed
     */
    final ButtonListener buttonListener = new ButtonListener();
        
    /**
     * Handle to the radio button that is currently selected
     */
    static MyRadioButton currentButton;
    
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

    /**
     * Constructor that creates all of the left pane's radiobuttons and JButtons, as well as, the right 
     * pane's scrollpane and JList.  It also creates all the listeners for these swing components so that
     * the simulation can be updated as needed.
     */
    public SimulationEditorPane(){
        setSize(splitPaneWidth, splitPaneHeight);
    	setLeftComponent(leftPane);
        setRightComponent(rightPane);
        setDividerLocation(0.4);
        leftPanelHeight = splitPaneHeight;
        leftPanelWidth = (int)(0.4*splitPaneWidth)-10;//10 accounts for width of slider bar
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 3;
        gbc.anchor = gbc.WEST;
        
        /**
         * This section sets up the JList of the right pane.  When an object from the list is selected,
         * its properties are displayed in the property sheet. 
         */
		rightPaneList.getSelectionModel().setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
		rightPaneList.addListSelectionListener(new MyListSelectionListener(){
		    public void valueChanged(javax.swing.event.ListSelectionEvent lse){
		        EditActions.setObject(rightPaneList.getSelectedValue());
		        setCurrentSelection(rightPaneList.getLeadSelectionIndex());
                if (added == false)
                    ViewActions.PROPERTYLIST.actionPerformed(new ActionEvent(this, 0, ""));
	            if (rightPaneList.getSelectedValue() != null){
                    wrapper = new Wrapper(rightPaneList.getSelectedValue(), rightPaneList.getSelectedValue().toString(), "simulate.gui." + rightPaneList.getSelectedValue().toString());
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
    	    }});// end JList instantiation and setup

        leftPanePanel.setLayout(gbl);
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
    }
    
    private void writeObject(ObjectOutputStream out) throws java.io.IOException {
        out.defaultWriteObject();
        out.writeObject(currentButton);
    }
    
    private void readObject(ObjectInputStream in) throws java.io.IOException, java.lang.ClassNotFoundException {
        in.defaultReadObject();
        currentButton = ((MyRadioButton)in.readObject());
        if (getTitle() == "Species")
            ViewActions.PROPERTYLIST.actionPerformed(new ActionEvent(this, 0, ""));
        wrapper = new Wrapper(new Blank(),"","");
        propertySheet.setTarget(wrapper);
    }    
    
    public void makeRadioButtons(Class[] className) {
		/**
		 * This section creates all of the left pane's radio buttons, adds listeners to these buttons,
		 * and adds them to the pane in a grid bag layout.
		 */
		leftPanelHeight = 0;
		for(int i=0; i < className.length; i++) {
            String name = className[i].getName();
            int idx;
            if (getTitle() == "Controller")
                idx = 8;
            else idx = 8 + getTitle().length();//name.indexOf("r");  //strip off simulate.P1 prefix
            name = name.substring(idx+1);            
            MyRadioButton button = new MyRadioButton(name,false,className[i]);
            button.addActionListener(new ButtonListener());
            gbl.setConstraints(button, gbc);
            simComponents.add(button);
            leftPanePanel.add(button);
            leftPanelHeight += radioButtonHeight;
            /* Use to make two columns out of the radiobuttons so that the panel isn't so tall
            if (className.length > 17) {
                leftPanelHeight -= 12;
                if (gbc.gridx == 0)
                    gbc.gridx = 3;
                else { 
                    gbc.gridx = 0;
                    gbc.gridy++;
                }
            }
            else*/ gbc.gridy++;
            
        }// end radio button addition
        if ((gbc.gridx == 3) && (className.length > 15)) gbc.gridy++;
        leftPanelHeight += 2*jButtonHeight;//Account for buttons on bottem of panel;
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
    }
    
    public final void checkSimFeasibility(){
        // Check if a sufficient number of components are added to allow a working simulation.
        // If so, enable the start button.
	    if (SimEditorTabMenu.allRemoveEnabled())
	        SimEditorTabMenu.setAllStart(true);
	    else SimEditorTabMenu.setAllStart(false);
    }
    
    public final void accountForNewSpecies(int i) {
        SpeciesPotentialLinkPane.setNumOSpecies(SpeciesPotentialLinkPane.getNumOSpecies() + i);
        SimEditorTabMenu.potential1Editor.update();
        SimEditorTabMenu.potential2Editor.update();
    }
    
    /**
     * This section creates the "Add" button which when pressed creates an instance of the class that
     * corresponds to the currently selected radio button of the left pane.  It then adds this object
     * to the JList of the rightpane and the simulation.instance object.
     */ 
    public final void addAddButton(){
        gbc.gridx = 0;
        gbc.gridwidth = 1;
        gbl.setConstraints(addToSim, gbc);
        if (title == "Phase")
            addToSim.setText("Add Phase");
        leftPanePanel.add(addToSim);
    }// end of add button
    
    /**
     * This section creates the start button that makes the internal frame that contains the 
     * simulation.instance object visible.  
     */
    public final void addStartButton(){
        gbc.gridx++;
	    start.setEnabled(false);
	    start.addActionListener(new MyActionListener(){
	        public void actionPerformed(ActionEvent e){
	            SimulateActions.getApplet().setVisible(true);
//                Simulation.instance.elementCoordinator.go();
                SimulateActions.getApplet().getContentPane().repaint();
                try { SimulateActions.getApplet().setSelected(true); }
                catch (java.beans.PropertyVetoException exc){} // attempt was vetoed
            }});
        gbl.setConstraints(start, gbc);
	    leftPanePanel.add(start);
    }// end of start button
    
	/**
     * This section creates a remove button which allows the user to remove components that are listed
	 * in the JList of the right pane.  If pressed, the currently selected object of the JList is 
	 * from both the JList and the simulation.instance object
	 */
    public final void addRemoveButton(){
	    gbc.gridx++;
	    gbl.setConstraints(remove, gbc);
	    leftPanePanel.add(remove);
    }// end of remove button
    
	/**
	 * This section creates the property sheet button which displays the properties of the currently
	 * selected object of the JList or null.
	 */
    public final void addPropertyButton(){
	    gbc.gridx = 0;
	    gbc.gridy++;
	    gbc.gridwidth = 3;
	    JButton propSheet = new JButton("Property Sheet");
	    propSheet.addActionListener(new MyActionListener(){
	        public void actionPerformed(ActionEvent e){
                if (added == false)
                    ViewActions.PROPERTYLIST.actionPerformed(new ActionEvent(this, 0, ""));
                if (rightPaneList.getSelectedValue() != null){
                    wrapper = new Wrapper(rightPaneList.getSelectedValue(), getTitle(), "simulate.gui." + getTitle()); 
                    propertySheet.setTarget(wrapper);
                }
                else {
                    wrapper = new Wrapper(FileActions.LOAD, "null", "null"); 
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
    }// end of property sheet button
    
    /**
     * Provides a handle to the single instance of the property sheet so that it can be updated with the
     * new component's properties.
     */
    public static void setPropertySheet(PropertySheet p){ propertySheet = p; }
    
    /**
     * Provides a handle to the wrapper object that is displayed in the propertySheet
     */
    public static void setWrapper(Wrapper w){ wrapper = w; }
    
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
     * Extension class of JRadioButton that allows button to know what simulation class they correspond to.
     */
    protected class MyRadioButton extends JRadioButton {
        Object cls;
        
        MyRadioButton(String label, boolean selected, Object cls) {
            super(label,selected);
            this.cls = cls;
        }// end of MyRadioButton constructor
    }// end of MyRadioButton class
    
    /**
     * Extension class of ActionListener that saves a handle to the radio button that is currently selected
     */
	protected class ButtonListener implements ActionListener, java.io.Serializable {
	    
	    public void actionPerformed(ActionEvent evt) {
	        currentButton = ((MyRadioButton)evt.getSource());
        }//end of actionPerformed
	}//end of ButtonListener
	
	public javax.swing.DefaultListModel getComponentList(){ return componentList; }
	
	private class MyListSelectionListener implements javax.swing.event.ListSelectionListener, java.io.Serializable {
	    public void valueChanged(javax.swing.event.ListSelectionEvent lse){}
	}
	
	public class MyInternalFrameAdapter extends javax.swing.event.InternalFrameAdapter implements java.io.Serializable {}

	public class MyActionListener implements java.awt.event.ActionListener, java.io.Serializable {
	    public void actionPerformed(java.awt.event.ActionEvent ae){}
	}
}//end of SimulationEditorPane class
