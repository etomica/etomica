/**
 * ControllerEditorPane
 *
 * This class extends the SimulationEditorPane and provides a point and click way of adding
 * controllers to a simulation.  It automatically locates all controllers that are present in the
 * system directory and displays them as JRadioButtons on a JPanel.  Once controllers are added,
 * they can be removed with the Remove button.  If enough components are present to start a 
 * simulation, the simulation can be instantiated by pushing the Start button.  The list of current
 * added controllers is displayed as a JList and if any of these added controllers are selected,
 * their properties are displayed in a PropertySheet.
 *
 * @author Bryan C. Mihalick
 * 10/10/00
 */
package etomica.gui;

import etomica.*;
import javax.swing.JButton;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class ControllerEditorPane extends SimulationEditorPane {
    private static int IDnumber = 0;
    
    ControllerEditorPane(SimulationEditor ed){
        super(ed);
        setTitle("Controller");
        
        // Calls introspection method from SimulationEditorPane to instantiate the related species
        // classes contained in the SimulateActions.controllerClasses static array.  These are listed
        // as a buttonGroup of radioButtons
        makeRadioButtons(SimulateActions.controllerClasses);
        remove.setEnabled(false);
		
        // The new actionListener will add an instance of the class that corresponds to the currently
        // selected radioButton to the simulationEditor.getSimulation() object as well as to the controllerEditorPane's
        // componentList.  
        addToSim.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    remove.setEnabled(true);                            // Enable 'Remove' button  

	                try {   // Try to make an instance of the selected class
	                    setComponent(((Class)currentButton.cls).newInstance());
	                    ((Controller)getComponent()).setName(((Class)currentButton.cls).getName().substring(8) + Integer.toString(IDnumber++));
	                    componentList.addElement(getComponent()); // Add new object to the componentList
                        simulationEditor.getSimulation().elementCoordinator.add((Simulation.Element)getComponent());
                        etomica.Device button = ((Controller)getComponent()).getButton();
                        if(button != null) simulationEditor.getSimulation().elementCoordinator.add(button);
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
                simulationEditor.getSimulation().unregister(((Controller)componentList.getElementAt(getCurrentSelection())));
	            componentList.remove(getCurrentSelection());
                wrapper = new Wrapper(FileActions.OPEN, "null", "null"); 
                propertySheet.setTarget(wrapper);

                if (componentList.getSize() == 0)
                    ((JButton)e.getSource()).setEnabled(false);
                
                simulationEditor.checkSimFeasibility();
	        }});
        addRemoveButton();  // Creates and adds the new JButton 'Remove'
        addPropertyButton();// Creates and adds the new JButton 'Property Sheet'
    }// end of ControllerEditorPane constructor
}// end of ControllerEditorPane class
        