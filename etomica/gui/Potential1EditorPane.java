package etomica.gui;

import etomica.Simulation;
import etomica.Potential1;
import etomica.Species;
import etomica.utility.HashMap2;
import java.awt.Color;
import javax.swing.JButton;
import javax.swing.JLabel;

//Java2 imports
//import java.util.LinkedList;

import etomica.utility.LinkedList;


public class Potential1EditorPane extends PotentialEditorPane {
    
    public HashMap2 potentialButtons = new HashMap2();
    
    Potential1EditorPane(SimulationEditor ed){
        super(ed);
        setTitle("Potential1");
    }
    
    public void update(){
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 1;
        gbc.anchor = gbc.WEST;
        leftPanelHeight = (speciesCount()+3)*jButtonHeight;
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
        leftPanePanel.removeAll();
        leftPanePanel.setLayout(gbl);
        JLabel label1 = new JLabel("   ");
        leftPanePanel.add(label1);
        gbc.gridx++;
        
        int potButtonsSpacer = 0;
        
        /**
         * P1 interactions are intramolecular and therefore don't need an upper triangular matrix
         * of buttons, unlike P2.  This accomplishes this single column of JButtons based on the 
         * numOSpecies.
         */
        currentButtons = new JButton[speciesCount()+1];
        potButtons = new JButton[speciesCount()+1];
 
        gbc.gridx = 0;
        gbc.gridy++;
        // add buttons and labels
        for (int i = 0; i < speciesCount(); i++){
            JLabel label2 = new JLabel("Species " + String.valueOf(i));
            gbl.setConstraints(label2,gbc);
            leftPanePanel.add(label2);
            gbc.gridx++;

            SpeciesPairButton button = new SpeciesPairButton(String.valueOf(i));
            potentialButtons.put(Integer.toString(0),Integer.toString(i),button);
            button.species1 = (Species)simulationEditor.getSimulation().speciesList().get(i);
            potButtons[i] = button;
            button.setBackground(Color.lightGray);
            button.addActionListener(new ButtonListener());
            gbl.setConstraints(button,gbc);
            leftPanePanel.add(button);
            gbc.gridx = 0;
            gbc.gridy++;
        }// end of button and label addition
        
        // As long as at least one species exists, the "Add," "Start," and "Remove" buttons are added
        if (speciesCount() != 0)
            addButtons();
        
//        linkButtons(); // Associate a potential to each button if one is present in the simulation
    }// end of update method
    
/*    private void linkButtons(){
        LinkedList list = simulationEditor().getSimulation().potential1List();

        for(int i = 0; i < list.size(); i++){
            Potential1 potential = (Potential1)list.get(i);
            SpeciesPairButton button = (SpeciesPairButton)potentialButtons.get(Integer.toString(0),Integer.toString(potential.getSpeciesIndex()));
            button.potential = potential.getClass();
            button.setEnabled(false);
            button.setBackground(Color.lightGray);
        }
    }
*/    
    public HashMap2 potentialButtons(){ return potentialButtons; }
}
