/**
 * PotentialEditorPane
 *
 * The PotentialEditorPane class is a splitpane that lists all the simulation components of a 
 * respective potential category (potential 1 or potential 2) on the leftside, and all of the added 
 * components from the leftside list on the rightside in a JList format.
 *
 * @author Bryan C. Mihalick
 * 8/16/00
 */

package etomica.gui;

import etomica.*;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.*;
import java.io.File;
import java.io.FilenameFilter;
import javax.swing.JButton;
import javax.swing.JRadioButton;
import javax.swing.event.InternalFrameEvent;
import java.util.jar.JarFile;
import java.util.zip.ZipEntry;

public abstract class PotentialEditorPane extends EditorPane {

    /**
     * Array of all current JButtons that are pressed 
     */
    public JButton[] currentButtons = new JButton[1];
    
    /**
     * Array of all potential JButtons on the left pane
     */
    public JButton[] potButtons = new JButton[1];
    
    /**
     * Number of buttons on left pane that are currently pressed
     */
    int buttonCount = 0;
    
    /**
     * Index of the potential class that the currently selected JButton in the left pane corresponds to
     */
    private int currentIndex = 0;
   
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
     * Determines if Potential Arrays need to be instantiated or not
     */
    public static boolean makePotArrays = true;
    
    /**
     * Makes it possible to determine which radio button was selected when the "add" button is pressed
     */
    final ButtonListener buttonListener = new ButtonListener();
    
    /**
     * Constructor that creates all of the left pane's radiobuttons and JButtons, as well as, the right 
     * pane's scrollpane and JList.  It also creates all the listeners for these swing components so that
     * the simulation can be updated as needed.
     */
    public PotentialEditorPane(SimulationEditor ed){
        super(ed);
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
		        Object obj = rightPaneList.getSelectedValue();
		        EditActions.setObject(obj);
		        
		        // See if customizer exists for the selected object.  If one does, enable customize
		        // selection on EditMenu, otherwise disable it.
		        try {
		            if (Introspector.getBeanInfo(obj.getClass()).getBeanDescriptor().getCustomizerClass() != null) {
		                EtomicaMenuBar.customizeItem.setEnabled(true);
		            }
		            else EtomicaMenuBar.customizeItem.setEnabled(false);
                }
                catch (IntrospectionException ie){}
		        setCurrentSelection(rightPaneList.getLeadSelectionIndex());
                for (int i = 0; i < leftPanePanel.getComponentCount(); i++){
                    if (leftPanePanel.getComponent(i) instanceof SpeciesPairButton){
                        if (((SpeciesPairButton)leftPanePanel.getComponent(i)).index == getCurrentSelection() && getCurrentSelection() != -1)
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
	                    Simulation.Element wrapperPot;
	                    if (rightPaneList.getSelectedValue() instanceof P2SimpleWrapper)
    	                    wrapperPot = ((P2SimpleWrapper)rightPaneList.getSelectedValue()).getOnlyPotential();
                        else wrapperPot = (Simulation.Element)rightPaneList.getSelectedValue();
                        propertySheet.setTarget(wrapperPot);
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
	                    if (((SpeciesPairButton)potButtons[i]).index == getCurrentSelection()){
	                        potButtons[i].setEnabled(true);
	                        potButtons[i].setBackground(Color.lightGray);
	                        ((SpeciesPairButton)potButtons[i]).index = -1;
	                        ((SpeciesPairButton)potButtons[i]).potential = null;
	                    }
	                }
	                currentIndex--;                         // Decreases the total number of current potentials
	                shiftIndices(getCurrentSelection());         // Updates indices of each button
	                if (getTitle() == "Potential1") {
	                        Simulation.instance.unregister((Potential1)componentList.getElementAt(getCurrentSelection()));
                    }
                    else {
                        Simulation.instance.unregister((Potential2)componentList.getElementAt(getCurrentSelection()));
	                }
	                componentList.remove(getCurrentSelection()); // Removes selected potential from the component list
                    propertySheet.setTarget(null);
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
                    Object obj = rightPaneList.getSelectedValue();
                    propertySheet.setTarget((Simulation.Element)obj);// Updates property sheet with the selected object's properties

		            // See if customizer exists for the selected object.  If one does, enable customize
		            // selection on EditMenu, otherwise disable it.
		            try {
		                if (Introspector.getBeanInfo(obj.getClass()).getBeanDescriptor().getCustomizerClass() != null) {
		                    EtomicaMenuBar.customizeItem.setEnabled(true);
		                }
		                else EtomicaMenuBar.customizeItem.setEnabled(false);
                    }
                    catch (IntrospectionException ie){}
                }
                else propertySheet.setTarget(null);
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
    
    public void setCurrentIndex(int cs) { currentIndex = cs; }
    public int getCurrentIndex() { return currentIndex; }
    
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
            setCurrentSelection(-1);
	        rightPaneList.clearSelection();
	        ((java.awt.Component)evt.getSource()).setBackground(Color.red);

            if (makePotArrays) {
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
	                potentialFrame.setVisible(true);
	                potentialFrame.setSelected(true);
	            }
	            catch(java.beans.PropertyVetoException pve){}
            }
        }// end of actionPerformed
	}// end of ButtonListener class
}// end of PotentialEditorPane class