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
 
package etomica.gui;

import etomica.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JRadioButton;
import javax.swing.event.InternalFrameEvent;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.File;
import java.io.FilenameFilter;
import java.beans.*;

public class SimulationEditorPane extends EditorPane {
    
    /**
     * Handle to the radio button that is currently selected
     */
    MyRadioButton currentButton;
    
    /**
     * Makes it possible to determine which radio button was selected when the "add" button is pressed
     */
    final ButtonListener buttonListener = new ButtonListener();
    
    private java.util.HashMap elementClasses = new java.util.HashMap(16);
    
    private int IDnumber = 0;
    
    private double dividerPlacement = 0.4;

    /**
     * Static array of all simulation components that extend species.class
     */
    public static Class[] speciesClasses;
    
    /**
     * Static array of all simulation components that extend integrator.class
     */
    public static Class[] integratorClasses;
    
    /**
     * Static array of all simulation components that extend phase.class
     */
    public static Class[] phaseClasses;
    
    /**
     * Static array of all simulation components that extend controller.class
     */
    public static Class[] controllerClasses;
    
    /**
     * Static array of all simulation components that extend display.class
     */
    public static Class[] displayClasses;
    
    /**
     * Static array of all simulation components that extend meter.class
     */
    public static Class[] meterClasses;
    
    /**
     * Static array of all simulation components that extend device.class
     */
    public static Class[] deviceClasses;

    /**
     * Constructor that creates all of the left pane's radiobuttons and JButtons, as well as, the right 
     * pane's scrollpane and JList.  It also creates all the listeners for these swing components so that
     * the simulation can be updated as needed.
     */
    public SimulationEditorPane(SimulationEditor ed, String t){
        super(ed);
        speciesClasses = IntrospectionArrays.speciesClasses;
        integratorClasses = IntrospectionArrays.integratorClasses;
        phaseClasses = IntrospectionArrays.phaseClasses;
        controllerClasses = IntrospectionArrays.controllerClasses;
        displayClasses = IntrospectionArrays.displayClasses;
        meterClasses = IntrospectionArrays.meterClasses;
        deviceClasses = IntrospectionArrays.deviceClasses;
        
        elementClasses.put("Species", speciesClasses);
        elementClasses.put("Integrator", integratorClasses);
        elementClasses.put("Phase", phaseClasses);
        elementClasses.put("Controller", controllerClasses);
        elementClasses.put("Display", displayClasses);
        elementClasses.put("Meter", meterClasses);
        elementClasses.put("Device", deviceClasses);
        
        setTitle(t);
        setSize(splitPaneWidth, splitPaneHeight);
    	setLeftComponent(leftPane);
        setRightComponent(rightPane);
	    if (getTitle() == "Phase") dividerPlacement = 0.5;
        else if (getTitle() == "Meter") dividerPlacement = 0.65;
	    setDividerLocation(dividerPlacement);
        leftPanelWidth = (int)(dividerPlacement*splitPaneWidth)-10;//10 accounts for width of slider bar
        leftPanelHeight = splitPaneHeight;
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 3;
        gbc.anchor = gbc.WEST;
        
        /**
         * This section sets up the JList of the right pane.  When an object from the list is selected,
         * its properties are displayed in the property sheet. 
         */
		rightPaneList.getSelectionModel().setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
		rightPaneList.addListSelectionListener(new MyListSelectionListener() {
		    public void valueChanged(javax.swing.event.ListSelectionEvent lse){
		        Object obj = rightPaneList.getSelectedValue();
		        EditActions.setObject(obj);
		        if(obj == null) return;
		        
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
                if (added == false)
                    ViewActions.PROPERTYLIST.actionPerformed(new ActionEvent(this, 0, ""));

                propertySheet.setTarget((Simulation.Element)obj);
                try {
                    propertySheet.setClosed(false);
                    propertySheet.setSelected(true);
                    propertySheet.setVisible(true);
                }
                catch(java.beans.PropertyVetoException pve){}
    	    }});// end of addListSelectionListener statement

	    propertySheet.addInternalFrameListener(new MyInternalFrameAdapter(){
	        public void internalFrameClosed( InternalFrameEvent ife ){
	            rightPaneList.clearSelection();
	        }});


        leftPanePanel.setLayout(gbl);
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
        /*if (getTitle() = "Meter"){
            leftPanePanel.setMinimumSize(new java.awt.Dimension(360, getHeight()));
	        leftPanePanel.setPreferredSize(new java.awt.Dimension(360, getHeight()));
        }*/
        
        // Calls introspection method from SimulationEditorPane to instantiate the related species
        // classes contained in the SimulateActions.<component name>Classes static array.  These are listed
        // as a buttonGroup of radioButtons
        Class[] classArray = (Class[])elementClasses.get(getTitle());
        if (getTitle() != "Phase") makeRadioButtons(classArray);
        remove.setEnabled(false);
		if (getTitle() == "Species") makeDefineMoleculeButton();
		
        // The new actionListener will add an instance of the class that corresponds to the currently
        // selected radioButton to the simulationEditor.getSimulation() object as well as to the SimulationEditorPane's
        // componentList and the javaWriter for the simulation.  
        addToSim.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    remove.setEnabled(true);                            // Enable 'Remove' button  

                    if ((getTitle() == "Phase") || (currentButton.cls != null)) {
	                    try {   // Try to make an instance of the selected class
	                        if (getTitle() == "Phase") {
	                            Phase phase = new Phase();
	                            setComponent(phase);
	                            ((Simulation.Element)getComponent()).setName(getTitle() + Integer.toString(IDnumber++));
	                        }
	                        else {
	                            setComponent(((Class)currentButton.cls).newInstance());
	                            ((Simulation.Element)getComponent()).setName(((Class)currentButton.cls).getName().substring(8) + Integer.toString(IDnumber++));
    	                    }
    	                    Simulation currentSimulation = simulationEditor.getSimulation();
    	                    Simulation.Element element = (Simulation.Element)getComponent();
	                        componentList.addElement(element); // Add new object to the componentList
                            currentSimulation.elementCoordinator.add(element);
                            ((JavaWriter)Etomica.javaWriters.get(currentSimulation)).add(element);
	                        if (getTitle() == "Controller"){
                                etomica.Device button = ((Controller)getComponent()).getButton();
                                if(button != null) simulationEditor.getSimulation().elementCoordinator.add(button);
	                        }
	                    }
	                    catch(InstantiationException exc) {}
	                    catch(IllegalAccessException exc) {}
	                }
                    else {
                        DefineMoleculeFrame newMolFrame = new DefineMoleculeFrame(simulationEditor);
                        Etomica.DesktopFrame.desktop.add(newMolFrame);
                        try { newMolFrame.setSelected(true); }
                        catch (java.beans.PropertyVetoException pve){}
                    }
                    
                    if (getTitle() == "Species") accountForNewSpecies();
                    // Check if the addition of the new species will complete the list of necessary
                    // components for a feasible simulation.  If so, it will enable the 'start' button.
                    simulationEditor.checkSimFeasibility();
             }});
        addAddButton();     // Creates and adds the new JButton 'Add'
        addStartButton();   // Creates and adds the new JButton 'Start'

        // The new actionListener will remove the object corresponding to the current selection of the
        // componentList from the simulationEditor.getSimulation() object.
	    remove.addActionListener(new MyActionListener(){
	        public void actionPerformed(ActionEvent e){
    	        Simulation currentSimulation = simulationEditor.getSimulation();
    	        Simulation.Element element = (Simulation.Element)componentList.getElementAt(getCurrentSelection());
                currentSimulation.unregister(element);
                ((JavaWriter)Etomica.javaWriters.get(currentSimulation)).remove(element);
	            componentList.remove(getCurrentSelection());
                propertySheet.setTarget(null);

                if (getTitle() == "Species") accountForNewSpecies();
                if (componentList.getSize() == 0)
                    ((JButton)e.getSource()).setEnabled(false);
                
                simulationEditor.checkSimFeasibility();
	        }});
        addRemoveButton();  // Creates and adds the new JButton 'Remove'
        addPropertyButton();// Creates and adds the new JButton 'Property Sheet'
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
        propertySheet.setTarget(null);
    }    
    
    public void makeRadioButtons(Class[] className) {
		/**
		 * This section creates all of the left pane's radio buttons, adds listeners to these buttons,
		 * and adds them to the pane in a grid bag layout.
		 */
		if (className == null){
		    leftPanelHeight = 0;
            leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	        leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
		    return;
		}
		 
		leftPanelHeight = 0;
		for(int i=0; i < className.length; i++) {
            String name = className[i].getName();
            int idx;
            if (getTitle() == "Controller")
                idx = 7;
            else idx = 7 + getTitle().length();//strip off etomica. and the <class name> prefix
            name = name.substring(idx+1);            
            MyRadioButton button = new MyRadioButton(name,false,className[i]);
        //new stuff
            Class[] interfaces = className[i].getInterfaces();
            for (int j = 0; j < interfaces.length; j++) {
                if (interfaces[j].toString().equals(String.valueOf("interface etomica.EtomicaElement"))){
                    etomica.EtomicaInfo info = null;
                    try {
                        java.lang.reflect.Method method = className[i].getMethod("getEtomicaInfo",null);
                        info = (etomica.EtomicaInfo)method.invoke(className[i], null);
  //                      info.setDescription(className[i].toString().substring(14));
                    }
                    catch(java.lang.SecurityException se){System.out.println("security exception");}
                    catch(java.lang.IllegalAccessException iae){System.out.println("illegal access exception");}
                    catch(java.lang.IllegalArgumentException ia){System.out.println("illegal argument exception");}
                    catch(java.lang.reflect.InvocationTargetException ite){System.out.println("invocation target exception");}
                    catch(java.lang.NoSuchMethodException nsme){
                        button.setToolTipText("No description available");
                        continue;
                    }
                    button.setToolTipText(info.getDescription());
                }
            }
        //end new stuff
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
    
    public void makeDefineMoleculeButton(){
		// Adds a RadioButton that allows the user to make their own Molecule with Atoms of different 
		// AtomType.
		MyRadioButton defineMolecule = new MyRadioButton("Define Molecule",false,null);
		defineMolecule.addActionListener(buttonListener);
		gbl.setConstraints(defineMolecule, gbc);
		simComponents.add(defineMolecule);
		leftPanePanel.add(defineMolecule);
		gbc.gridy++;
        leftPanelHeight += radioButtonHeight;//Account for extra radioButton on leftPanePanel;
    }
    
    public final void accountForNewSpecies() {
        simulationEditor.potential1Editor.update();
        simulationEditor.potential2Editor.update();
    }
    
    /**
     * This section creates the "Add" button which when pressed creates an instance of the class that
     * corresponds to the currently selected radio button of the left pane.  It then adds this object
     * to the JList of the rightpane and the getSimulation() object.
     */ 
    public final void addAddButton(){
        gbc.gridx = 0;
        gbc.gridwidth = 1;
        gbl.setConstraints(addToSim, gbc);
        if (getTitle() == "Phase")
            addToSim.setText("Add Phase");
        leftPanePanel.add(addToSim);
    }// end of add button
    
    /**
     * This section creates the start button that makes the internal frame that contains the 
     * simulation object visible.  
     */
    public final void addStartButton(){
        gbc.gridx++;
	    start.setEnabled(false);
	    start.addActionListener(ViewActions.SIMULATION);/*new MyActionListener(){
	        public void actionPerformed(ActionEvent e){
	            SimulationFrame frame = (SimulationFrame)Etomica.simulationFrames.get(simulationEditor.getSimulation());
	            frame.setVisible(true);
                try { frame.setSelected(true); }
                catch (java.beans.PropertyVetoException exc){} // attempt was vetoed
 //               frame.getContentPane().repaint();
            }});*/
        gbl.setConstraints(start, gbc);
	    leftPanePanel.add(start);
    }// end of start button
    
	/**
     * This section creates a remove button which allows the user to remove components that are listed
	 * in the JList of the right pane.  If pressed, the currently selected object of the JList is 
	 * from both the JList and the simulation object
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
                    Object obj = rightPaneList.getSelectedValue();
                    propertySheet.setTarget((Simulation.Element)obj);

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
                else { propertySheet.setTarget(null); }
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
}//end of SimulationEditorPane class
