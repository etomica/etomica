/**
 * Potential Frame
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
import java.beans.PropertyVetoException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JInternalFrame;

public class AtomPotentialFrame extends PotentialFrame {
    /**
     * Constructor that takes the potential type (P1 or P2) and the removeButton from the Potential tab
     * of the SimEditorTabMenu that called this constructor.  It determines the necessary size based on
     * the type of potentials being used, sets the title, and adds the JInternalFrame to the JDesktopPane
     */
    public AtomPotentialFrame(){
        setBounds(515,60,250,420);
        setTitle("Atom Potentials");  
        panelWidth = 200;
        potentialPanel.setSize(250, 420);    
    
        makeRadioButtons(potentialClasses,17); 
            
        addToSim.addActionListener(new MyActionListener(){
            public void actionPerformed(ActionEvent e){
                int arrayIndex1 = ((DefineAtomPotentialFrame.SpeciesPairButton)DefineAtomPotentialFrame.currentButton).speciesIndex1;
                int arrayIndex2 = ((DefineAtomPotentialFrame.SpeciesPairButton)DefineAtomPotentialFrame.currentButton).speciesIndex2;
                        
                try {
                    atomPairPotArray[arrayIndex1][arrayIndex2] = ((Class)currentButton.cls);
   	            }
	            catch(java.lang.ArrayIndexOutOfBoundsException exc) {atomPairPotArray[arrayIndex1][arrayIndex1 + 1] = /*(Potential)*/((Class)currentButton.cls); }//.newInstance();}
                        
	            try {
                    ((JInternalFrame)Etomica.DesktopFrame.desktop.getComponent(0)).setClosed(true);
                }
                catch(PropertyVetoException pve){}

	            addToSim.removeActionListener(this);
            }});
        
        addAddButton();
        
        cancel.addActionListener(new MyActionListener(){
                public void actionPerformed(ActionEvent e){
                    for (int i = 0; getPotentialEditor().currentButtons[i] != null; i++){
                        getPotentialEditor().currentButtons[i].setBackground(Color.lightGray);
                        getPotentialEditor().currentButtons[i] = null;
                    }
                    getPotentialEditor().buttonCount = 0;
	                    
	                try {
                        ((JInternalFrame)Etomica.DesktopFrame.desktop.getComponent(0)).setClosed(true);
                    }
                    catch(PropertyVetoException pve){}                        
                }});
        addCancelButton();
        panelHeight += jButtonHeight;// Accounts for height of the jButtons that are added
        potentialPanel.setMinimumSize(new Dimension(panelWidth, panelHeight));
        potentialPanel.setPreferredSize(new Dimension(panelWidth, panelHeight));
        setSize(panelWidth+20,panelHeight+45);
    }// end of PotentialFramePanel
}// end of PotentialFrame class