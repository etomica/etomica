/**
 * PhaseEditorPane
 *
 * This class extends the SimulationEditorPane and provides a point and click way of adding
 * phases to a simulation.  It automatically locates all phases that are present in the
 * system directory and displays them as JRadioButtons on a JPanel.  Once a phase is added,
 * it can be removed with the Remove button.  If enough components are present to start a 
 * simulation, the simulation can be instantiated by pushing the Start button.  The list of current
 * added phases is displayed as a JList and if any of these added phases are selected,
 * their properties are displayed in a PropertySheet.
 *
 * @author Bryan C. Mihalick
 * 9/29/00
 */
 
package simulate.gui;

import simulate.*;
import javax.swing.JButton;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class PhaseEditorPane extends SimulationEditorPane {
    private static int IDnumber = 0;
    
    PhaseEditorPane(SimulationEditor ed){
        super(ed);
        setTitle("Phase");
        setDividerLocation(0.5);
        leftPanelWidth = (int)(0.5*splitPaneWidth)-10;
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
        
        // Calls introspection method from SimulationEditorPane to instantiate the related species
        // classes contained in the SimulateActions.phaseClasses static array.  These are listed
        // as a buttonGroup of radioButtons
        makeRadioButtons(SimulateActions.phaseClasses);
        remove.setEnabled(false);
		
        // This section makes sure that a phase is always added at the beginning by default since no
		// simulation can run without one.
		((javax.swing.AbstractButton)leftPanePanel.getComponent(0)).doClick();
        remove.setEnabled(true);
	    try {
	        setComponent(((Class)currentButton.cls).newInstance());
	        ((Phase)getComponent()).setName(((Class)currentButton.cls).getName().substring(9) + Integer.toString(IDnumber++));
	    }
	    catch(InstantiationException exc) {}
	    catch(IllegalAccessException exc) {}
        simulationEditor.getSimulation().elementCoordinator.add((Simulation.Element)getComponent());
        componentList.addElement(getComponent());
	    if (simulationEditor.potential2Editor.remove.isEnabled() && simulationEditor.speciesEditor.remove.isEnabled() && 
	        simulationEditor.integratorEditor.remove.isEnabled() && simulationEditor.phaseEditor.remove.isEnabled() && 
	        simulationEditor.controllerEditor.remove.isEnabled() && simulationEditor.displayEditor.remove.isEnabled()){
	            simulationEditor.setAllStart(true);
	    }
		leftPanePanel.getComponent(0).setVisible(false);
		
        // The new actionListener will add an instance of the class that corresponds to the currently
        // selected radioButton to the simulationEditor.getSimulation() object as well as to the phaseEditorPane's
        // componentList.  
        addToSim.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    EtomicaMenuBar.selectSpaceItem.setEnabled(false);   // Disable 'Select Space' menuItem
                    remove.setEnabled(true);                            // Enable 'Remove' button  

	                //try {   // Try to make an instance of the selected class
	                    setComponent(new Phase());//((Class)currentButton.cls).newInstance());
	                    ((Phase)getComponent()).setName("Phase" + Integer.toString(IDnumber++));
	                    componentList.addElement(getComponent()); // Add new object to the componentList
                        simulationEditor.getSimulation().elementCoordinator.add((Simulation.Element)getComponent());
	                //}
	                //catch(InstantiationException exc) {}
	                //catch(IllegalAccessException exc) {}

                    // Check if the addition of the new species will complete the list of necessary
                    // components for a feasible simulation.  If so, it will enable the 'start' button.
                    checkSimFeasibility();
             }});
        addAddButton();     // Creates and adds the new JButton 'Add'
        addStartButton();   // Creates and adds the new JButton 'Start'

        // The new actionListener will remove the object corresponding to the current selection of the
        // componentList from the simulationEditor.getSimulation() object.
	    remove.addActionListener(new MyActionListener(){
	        public void actionPerformed(ActionEvent e){
                simulationEditor.getSimulation().unregister(((Phase)componentList.getElementAt(getCurrentSelection())));
	            componentList.remove(getCurrentSelection());
                wrapper = new Wrapper(FileActions.OPEN, "null", "null"); 
                propertySheet.setTarget(wrapper);

                if (componentList.getSize() == 0)
                    ((JButton)e.getSource()).setEnabled(false);
                
                checkSimFeasibility();
	        }});
        addRemoveButton();  // Creates and adds the new JButton 'Remove'
        addPropertyButton();// Creates and adds the new JButton 'Property Sheet'
    }// end of PhaseEditorPane constructor
}// end of PhaseEditorPane class
        