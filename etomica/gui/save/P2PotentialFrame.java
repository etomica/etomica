/**
 * P2Potential Frame
 *
 * This class makes the JInternalFrame that lists all the related P1 or P2 potentials from the simulation
 * source files.  It uses a buttongroup of JRadiobuttons to determine which potential is picked.  Once
 * picked, the selected potential is instantiated, added to the simulation.instance object, and added to
 * the component list of the "Potential" tab of the SimEditorTabMenu.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */

package simulate.gui;

import simulate.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JInternalFrame;

public class P2PotentialFrame extends PotentialFrame {
    private static int IDnumber = 0;
    private SimulationEditor simulationEditor;
    
    /**
     * Constructor that takes the potential type (P1 or P2) and the removeButton from the Potential tab
     * of the simulationEditor that called this constructor.  It determines the necessary size based on
     * the type of potentials being used, sets the title, and adds the JInternalFrame to the JDesktopPane
     */
    public P2PotentialFrame(SimulationEditor ed){
        simulationEditor = ed;
        setBounds(515,60,200,345);
        setTitle("P2 Potentials");  
        panelWidth = 200;
        potentialPanel.setSize(200, 420);

        makeRadioButtons(potential2Classes,17);
        //makeRadioButtons(potentialClasses,17);
        
        // Adds the RadioButton that allows you to make your own special Potentials
		MyRadioButton defineMolecule = new MyRadioButton("Define Potential", false, null);
		defineMolecule.addActionListener(buttonListener);
		gbl.setConstraints(defineMolecule, gbc);
		simComponents.add(defineMolecule);
		potentialPanel.add(defineMolecule);
		gbc.gridy++;
        panelHeight += radioButtonHeight;// Accounts for height of this radioButton
		
        addToSim.addActionListener(new MyActionListener(){
            public void actionPerformed(ActionEvent e){
                EtomicaMenuBar.selectSpaceItem.setEnabled(false);
                simulationEditor.getParent().repaint();
                        
                for (int i = 0; potentialEditor.currentButtons[i] != null; i++){
                    ((SpeciesPotentialLinkPane.SpeciesPairButton)potentialEditor.currentButtons[i]).potential = ((Class)currentButton.cls);
                    potentialEditor.currentButtons[i].setEnabled(false);
                    potentialEditor.currentButtons[i].setBackground(Color.lightGray);
                }
                        
                potentialEditor.setCurrentIndex(potentialEditor.getCurrentIndex() + 1);
                potentialEditor.buttonCount = 0;
                removeButton.setEnabled(true);

	            try {
                    ((JInternalFrame)Etomica.DesktopFrame.desktop.getComponent(0)).setClosed(true);
                }
                catch(java.beans.PropertyVetoException pve){}
	                    
	            if (instantiate == true){
                    if (currentButton.cls != null){
	                    try {
	                        if (currentButton.cls.toString().startsWith("class simulate.Potential"))
	                            component = new P2SimpleWrapper((Potential)((Class)currentButton.cls).newInstance());
	                        else component = ((Class)currentButton.cls).newInstance();

                            ((Potential2)component).setSpecies1Index(((SpeciesPotentialLinkPane.SpeciesPairButton)potentialEditor.currentButtons[0]).speciesIndex1);
                            ((Potential2)component).setSpecies2Index(((SpeciesPotentialLinkPane.SpeciesPairButton)potentialEditor.currentButtons[0]).speciesIndex2);
	                        ((Potential2)component).setName(((Class)currentButton.cls).getName().substring(9) + Integer.toString(IDnumber++));
                             Simulation.instance.elementCoordinator.add((Simulation.Element)component);
	                    }
	                    catch(InstantiationException exc) {System.out.println(e.toString()); System.exit(1);}
	                    catch(IllegalAccessException exc) {System.out.println(e.toString()); System.exit(1);}
  
                        for (int i = 0; potentialEditor.currentButtons[i] != null; i++){
                            Simulation.instance.potential2[((SpeciesPotentialLinkPane.SpeciesPairButton)potentialEditor.currentButtons[i]).speciesIndex1][((SpeciesPotentialLinkPane.SpeciesPairButton)potentialEditor.currentButtons[i]).speciesIndex2] = (Potential2)component;
                            Simulation.instance.potential2[((SpeciesPotentialLinkPane.SpeciesPairButton)potentialEditor.currentButtons[i]).speciesIndex2][((SpeciesPotentialLinkPane.SpeciesPairButton)potentialEditor.currentButtons[i]).speciesIndex1] = (Potential2)component;
                            potentialEditor.currentButtons[i] = null;
                        }
                        potentialEditor.componentList.addElement(component); 
                    }
                    else {
                        DefineAtomPotentialFrame newPotFrame = new DefineAtomPotentialFrame();
                        Etomica.DesktopFrame.desktop.add(newPotFrame);
                        try { newPotFrame.setSelected(true); }
                        catch (java.beans.PropertyVetoException pve){}
                    }
	            }
	                    
                for (int i = 0; potentialEditor.currentButtons[i] != null; i++){
	                potentialEditor.currentButtons[i] = null;
	            }
	                    

	            if (simulationEditor.allRemoveEnabled())
	                    simulationEditor.setAllStart(true);

	            addToSim.removeActionListener(this);
            }});
        addAddButton();
        
        cancel.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    for (int i = 0; potentialEditor.currentButtons[i] != null; i++){
                        potentialEditor.currentButtons[i].setBackground(Color.lightGray);
                        potentialEditor.currentButtons[i] = null;
                    }
                    potentialEditor.buttonCount = 0;
    	                    
	                try {
                        ((JInternalFrame)Etomica.DesktopFrame.desktop.getComponent(0)).setClosed(true);
                    }
                    catch(java.beans.PropertyVetoException pve){}                        
                }});
        addCancelButton();
        panelHeight += jButtonHeight;// Accounts for height of the jButtons that are added
        potentialPanel.setMinimumSize(new Dimension(panelWidth, panelHeight));
        potentialPanel.setPreferredSize(new Dimension(panelWidth, panelHeight));
        setSize(panelWidth+20,panelHeight+45);
    }// end of P2PotentialFrame constructor
}// end of P2PotentialFrame class