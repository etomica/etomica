/**
 * SpeciesEditorPane
 *
 * An extension of the SimulationEditorPane.  The SpeciesEditorPane allows users to add and remove Species
 * objects from the Simulation.instance object.  This object also provides a 'define Molecule' radioButton
 * that when selected gives users the flexibility to create their own molecules with their own chosen
 * atoms and configurations.
 *
 * @author Bryan C. Mihalick
 * 9/29/00
 */

package simulate.gui;

import simulate.*;
import javax.swing.JButton;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class SpeciesEditorPane extends SimulationEditorPane {
    private static int IDnumber = 0;
    
    SpeciesEditorPane(){
        super();
        setTitle("Species");
        
        // Calls introspection method from SimulationEditorPane to instantiate the related species
        // classes contained in the SimulateActions.speciesClasses static array.  These are listed
        // as a buttonGroup of radioButtons
        makeRadioButtons(SimulateActions.speciesClasses);
        remove.setEnabled(false);
		
		// Adds a RadioButton that allows the user to make their own Molecule with Atoms of different 
		// AtomType.
		MyRadioButton defineMolecule = new MyRadioButton("Define Molecule",false,null);
		defineMolecule.addActionListener(buttonListener);
		gbl.setConstraints(defineMolecule, gbc);
		simComponents.add(defineMolecule);
		leftPanePanel.add(defineMolecule);
		gbc.gridy++;
        leftPanelHeight += radioButtonHeight;//Account for extra radioButton on leftPanePanel;
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));

        // The new actionListener will add an instance of the class that corresponds to the currently
        // selected radioButton to the Simulation.instance object as well as to the speciesEditorPane's
        // componentList.  
        addToSim.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    EtomicaMenuBar.selectSpaceItem.setEnabled(false);   // Disable 'Select Space' menuItem
                    remove.setEnabled(true);                            // Enable 'Remove' button  
                    
                    // If currentButton.cls == null, that means the 'Define Molecule' radioButton was
                    // selected.  Therefore, the DefineMoleculeFrame must be opened.  Otherwise, add
                    // an instance of class corresponding to the currently selected radioButton.
                    if (currentButton.cls != null){
	                    try {   // Try to make an instance of the selected class
	                        setComponent(((Class)currentButton.cls).newInstance());
	                        ((Species)getComponent()).setName(((Class)currentButton.cls).getName().substring(9) + Integer.toString(IDnumber++));
	                        componentList.addElement(getComponent()); // Add new object to the componentList
                            Simulation.instance.elementCoordinator.add((Simulation.Element)getComponent());
	                    }
	                    catch(InstantiationException exc) {}
	                    catch(IllegalAccessException exc) {}
                    }
                    else {
                        SimEditorTabMenu.potential1Editor.setEnabled(false);
                        SimEditorTabMenu.potential2Editor.setEnabled(false);
                        SimEditorTabMenu.integratorEditor.setEnabled(false);
                        SimEditorTabMenu.phaseEditor.setEnabled(false);
                        SimEditorTabMenu.controllerEditor.setEnabled(false);
                        SimEditorTabMenu.displayEditor.setEnabled(false);
                        SimEditorTabMenu.deviceEditor.setEnabled(false);
                        SimEditorTabMenu.meterEditor.setEnabled(false);
                        DefineMoleculeFrame newMolFrame = new DefineMoleculeFrame();
                        Etomica.DesktopFrame.desktop.add(newMolFrame);
                        try { newMolFrame.setSelected(true); }
                        catch (java.beans.PropertyVetoException pve){}
                    }
                    // Update total number of species as well as the SpeciesPotentialLinkPanes to 
                    // account for these new species
                    accountForNewSpecies(1);
                    
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
                Simulation.instance.unregister(((Species)componentList.getElementAt(getCurrentSelection())));
	            componentList.remove(getCurrentSelection());
                wrapper = new Wrapper(FileActions.LOAD, "null", "null"); 
                propertySheet.setTarget(wrapper);

                // Update total number of species as well as the SpeciesPotentialLinkPanes to 
                // account for these new species
                accountForNewSpecies(-1);
                
                // Disable the 'remove' button if no other components exist in the componentList.
                if (componentList.getSize() == 0)
                    ((JButton)e.getSource()).setEnabled(false);
                
                // Check if the deletion of the species will make it impossible for a simulation to
                // exist.  If so, it will disable the 'start' button.
                checkSimFeasibility();
	        }});
        addRemoveButton();  // Creates and adds the new JButton 'Remove'
        addPropertyButton();// Creates and adds the new JButton 'Property Sheet'
    }// end of SpeciesEditorPane constructor
}// end of SpeciesEditorPane class
        