/**
 * MeterEditorPane
 *
 * This class extends the SimulationEditorPane and provides a point and click way of adding
 * meters to a simulation.  It automatically locates all meters that are present in the
 * system directory and displays them as JRadioButtons on a JPanel.  Once a meter is added,
 * it can be removed with the Remove button.  If enough components are present to start a 
 * simulation, the simulation can be instantiated by pushing the Start button.  The list of current
 * added meters is displayed as a JList and if any of these added meters are selected,
 * their properties are displayed in a PropertySheet.
 *
 * @author Bryan C. Mihalick
 * 9/29/00
 */
 
package etomica.gui;

import etomica.*;
import javax.swing.JButton;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class MeterEditorPane extends SimulationEditorPane {
    private static int IDnumber = 0;
    
    MeterEditorPane(SimulationEditor ed){
        super(ed);
        setTitle("Meter");
        setDividerLocation(0.65);
        leftPanePanel.setMinimumSize(new java.awt.Dimension(360, getHeight()));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(360, getHeight()));	    
	
        // Calls introspection method from SimulationEditorPane to instantiate the related species
        // classes contained in the SimulateActions.meterClasses static array.  These are listed
        // as a buttonGroup of radioButtons
        makeRadioButtons(SimulateActions.meterClasses);
        remove.setEnabled(false);
		
        // The new actionListener will add an instance of the class that corresponds to the currently
        // selected radioButton to the simulationEditor.getSimulation() object as well as to the meterEditorPane's
        // componentList.  
        addToSim.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    remove.setEnabled(true);                            // Enable 'Remove' button  

	                try {   // Try to make an instance of the selected class
	                    setComponent(((Class)currentButton.cls).newInstance());
	                    ((MeterAbstract)getComponent()).setName(((Class)currentButton.cls).getName().substring(8) + Integer.toString(IDnumber++));
	                    componentList.addElement(getComponent()); // Add new object to the componentList
                        simulationEditor.getSimulation().elementCoordinator.add((Simulation.Element)getComponent());
	                }
	                catch(InstantiationException exc) {}
	                catch(IllegalAccessException exc) {}

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
                simulationEditor.getSimulation().unregister(((Meter)componentList.getElementAt(getCurrentSelection())));
	            componentList.remove(getCurrentSelection());
                wrapper = new Wrapper(FileActions.OPEN, "null", "null"); 
                propertySheet.setTarget(wrapper);

                if (componentList.getSize() == 0)
                    ((JButton)e.getSource()).setEnabled(false);
                
                simulationEditor.checkSimFeasibility();
	        }});
        addRemoveButton();  // Creates and adds the new JButton 'Remove'
        addPropertyButton();// Creates and adds the new JButton 'Property Sheet'
    }// end of MeterEditorPane constructor
}// end of MeterEditorPane class
        