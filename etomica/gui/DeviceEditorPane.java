/**
 * DeviceEditorPane
 *
 * This class extends the SimulationEditorPane and provides a point and click way of adding
 * devices to a simulation.  It automatically locates all devices that are present in the
 * system directory and displays them as JRadioButtons on a JPanel.  Once an device is added,
 * it can be removed with the Remove button.  If enough components are present to start a 
 * simulation, the simulation can be instantiated by pushing the Start button.  The list of current
 * added devices is displayed as a JList and if any of these added devices are selected,
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

public class DeviceEditorPane extends SimulationEditorPane {
    private static int IDnumber = 0;
    
    DeviceEditorPane(){
        super();
        setTitle("Device");
        
        // Calls introspection method from SimulationEditorPane to instantiate the related species
        // classes contained in the SimulateActions.deviceClasses static array.  These are listed
        // as a buttonGroup of radioButtons
        makeRadioButtons(SimulateActions.deviceClasses);
        remove.setEnabled(false);
		
        // The new actionListener will add an instance of the class that corresponds to the currently
        // selected radioButton to the Simulation.instance object as well as to the deviceEditorPane's
        // componentList.  
        addToSim.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    EtomicaMenuBar.selectSpaceItem.setEnabled(false);   // Disable 'Select Space' menuItem
                    remove.setEnabled(true);                            // Enable 'Remove' button  

	                try {   // Try to make an instance of the selected class
	                    setComponent(((Class)currentButton.cls).newInstance());
	                    ((Device)getComponent()).setName(((Class)currentButton.cls).getName().substring(9) + Integer.toString(IDnumber++));
	                    componentList.addElement(getComponent()); // Add new object to the componentList
                        Simulation.instance.elementCoordinator.add((Simulation.Element)getComponent());
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
        // componentList from the Simulation.instance object.
	    remove.addActionListener(new MyActionListener(){
	        public void actionPerformed(ActionEvent e){
                Simulation.instance.unregister(((Device)componentList.getElementAt(getCurrentSelection())));
	            componentList.remove(getCurrentSelection());
                wrapper = new Wrapper(FileActions.LOAD, "null", "null"); 
                propertySheet.setTarget(wrapper);

                if (componentList.getSize() == 0)
                    ((JButton)e.getSource()).setEnabled(false);
                
                checkSimFeasibility();
	        }});
        addRemoveButton();  // Creates and adds the new JButton 'Remove'
        addPropertyButton();// Creates and adds the new JButton 'Property Sheet'
    }// end of DeviceEditorPane constructor
}// end of DeviceEditorPane class
        