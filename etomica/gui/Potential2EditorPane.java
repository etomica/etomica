package etomica.gui;

import etomica.Simulation;
import etomica.Potential2;
import etomica.utility.HashMap2;
import java.awt.Color;
import javax.swing.JButton;
import javax.swing.JLabel;
//Java2 imports
//import java.util.LinkedList;

import etomica.utility.LinkedList;

public class Potential2EditorPane extends PotentialEditorPane {
    HashMap2 potentialButtons = new HashMap2();
    
    Potential2EditorPane(SimulationEditor ed){ 
        super(ed);
        setTitle("Potential2");
    }
    
    public void update(){
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 1;
        gbc.anchor = gbc.WEST;
        leftPanelHeight = 0;
        leftPanePanel.removeAll();
        leftPanePanel.setLayout(gbl);
        JLabel label1 = new JLabel("   ");
        leftPanePanel.add(label1);
        gbc.gridx++;
        
        int potButtonsSpacer = 0;
        
        /**
         * For P2 we have a different story since these interactions are intermolecular.  This section
         * accomplishes this.
         */
        currentButtons = new JButton[speciesCount()*(speciesCount()+1)/2+1];
        potButtons = new JButton[speciesCount()*(speciesCount()+1)/2+1];
        for (int i = 0; i < speciesCount(); i++){
            JLabel label2 = new JLabel("Species " + String.valueOf(i));
            gbl.setConstraints(label2,gbc);
            leftPanePanel.add(label2);
            gbc.gridx++;
        }
        leftPanelHeight += radioButtonHeight;//Account for row of labels;
        
        gbc.gridx = 0;
        gbc.gridy++;
        for (int i = 0; i < speciesCount(); i++){
            leftPanelHeight += jButtonHeight;// Account for each new row of JButtons
            JLabel label3 = new JLabel("Species " + String.valueOf(i));
            gbl.setConstraints(label3,gbc);
            leftPanePanel.add(label3);
            gbc.gridx++;
            // add labels
            for (int j = speciesCount()-i; j < speciesCount(); j++){
                JLabel label2 = new JLabel("   ");
                gbl.setConstraints(label2,gbc);
                leftPanePanel.add(label2);
                gbc.gridx++;
            }
            // add buttons with listeners and change color
            for (int k = 0; k < speciesCount()-i; k++){
                SpeciesPairButton button = new SpeciesPairButton(String.valueOf(i) + "," + String.valueOf(k+i));
                potentialButtons.put(Integer.toString(i),Integer.toString(k+i),button);
                button.speciesIndex1 = i;
                button.speciesIndex2 = k+i;
                potButtons[k + potButtonsSpacer] = button;
                button.setBackground(Color.lightGray);
                button.addActionListener(buttonListener);
                gbl.setConstraints(button,gbc);
                leftPanePanel.add(button);
                gbc.gridx++;
            }
            potButtonsSpacer += speciesCount()-i;
            gbc.gridx = 0;
            gbc.gridy++;
        }// end of button and label addition
           
        // As long as at least one species exists, the "Add," "Start," and "Remove" buttons are added
        if (speciesCount() != 0) 
            addButtons();
        
        leftPanelHeight += 2*jButtonHeight;// Accounts for "Add," "Start," "Remove," and "Property Sheet" Buttons
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
    
        linkButtons();
    }
    
    private void linkButtons(){
        LinkedList list = simulationEditor().getSimulation().potential2List();

        for(int i = 0; i < list.size(); i++){
            Potential2 potential = (Potential2)list.get(i);
            SpeciesPairButton button = (SpeciesPairButton)potentialButtons.get(Integer.toString(potential.getSpecies1Index()),Integer.toString(potential.getSpecies2Index()));
            button.potential = potential.getClass();
            button.setEnabled(false);
            button.setBackground(Color.lightGray);
        }
    }

    public HashMap2 potentialButtons(){ return potentialButtons; }
}