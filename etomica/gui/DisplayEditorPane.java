/**
 * DisplayEditorPane
 *
 * This class extends the SimulationEditorPane and provides a point and click way of adding
 * displays to a simulation.  It automatically locates all displays that are present in the
 * system directory and displays them as JRadioButtons on a JPanel.  Once an display is added,
 * it can be removed with the Remove button.  If enough components are present to start a 
 * simulation, the simulation can be instantiated by pushing the Start button.  The list of current
 * added displays is displayed as a JList and if any of these added displays are selected,
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

public class DisplayEditorPane extends SimulationEditorPane {
    private static int IDnumber = 0;
    
    DisplayEditorPane(SimulationEditor ed){
        super(ed);
        setTitle("Display");
        
        // Calls introspection method from SimulationEditorPane to instantiate the related species
        // classes contained in the SimulateActions.displayClasses static array.  These are listed
        // as a buttonGroup of radioButtons
        makeRadioButtons(SimulateActions.displayClasses);
        remove.setEnabled(false);
		
        // The new actionListener will add an instance of the class that corresponds to the currently
        // selected radioButton to the simulationEditor.getSimulation() object as well as to the displayEditorPane's
        // componentList.  
        addToSim.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    remove.setEnabled(true);                            // Enable 'Remove' button  

	                try {   // Try to make an instance of the selected class
	                    setComponent(((Class)currentButton.cls).newInstance());
	                    ((Display)getComponent()).setName(((Class)currentButton.cls).getName().substring(9) + Integer.toString(IDnumber++));
	                    componentList.addElement(getComponent()); // Add new object to the componentList
                        simulationEditor.getSimulation().elementCoordinator.add((Simulation.Element)getComponent());
	                }
	                catch(InstantiationException exc) {}
	                catch(IllegalAccessException exc) {}

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
                simulationEditor.getSimulation().unregister(((Display)componentList.getElementAt(getCurrentSelection())));
	            componentList.remove(getCurrentSelection());
                wrapper = new Wrapper(FileActions.OPEN, "null", "null"); 
                propertySheet.setTarget(wrapper);

                if (componentList.getSize() == 0)
                    ((JButton)e.getSource()).setEnabled(false);
                
                checkSimFeasibility();
	        }});
        addRemoveButton();  // Creates and adds the new JButton 'Remove'
        addPropertyButton();// Creates and adds the new JButton 'Property Sheet'
    }// end of DisplayEditorPane constructor
}// end of DisplayEditorPane class
        