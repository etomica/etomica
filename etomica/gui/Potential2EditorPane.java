package simulate.gui;

import java.awt.Color;
import javax.swing.JButton;
import javax.swing.JLabel;

public class Potential2EditorPane extends SpeciesPotentialLinkPane {
    
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
        currentButtons = new JButton[getNumOSpecies()*(getNumOSpecies()+1)/2+1];
        potButtons = new JButton[getNumOSpecies()*(getNumOSpecies()+1)/2+1];
        for (int i = 0; i < getNumOSpecies(); i++){
            JLabel label2 = new JLabel("Species " + String.valueOf(i));
            gbl.setConstraints(label2,gbc);
            leftPanePanel.add(label2);
            gbc.gridx++;
        }
        leftPanelHeight += radioButtonHeight;//Account for row of labels;
        
        gbc.gridx = 0;
        gbc.gridy++;
        for (int i = 0; i < getNumOSpecies(); i++){
            leftPanelHeight += jButtonHeight;// Account for each new row of JButtons
            JLabel label3 = new JLabel("Species " + String.valueOf(i));
            gbl.setConstraints(label3,gbc);
            leftPanePanel.add(label3);
            gbc.gridx++;
            // add labels
            for (int j = getNumOSpecies()-i; j < getNumOSpecies(); j++){
                JLabel label2 = new JLabel("   ");
                gbl.setConstraints(label2,gbc);
                leftPanePanel.add(label2);
                gbc.gridx++;
            }
            // add buttons with listeners and change color
            for (int k = 0; k < getNumOSpecies()-i; k++){
                SpeciesPairButton button = new SpeciesPairButton(String.valueOf(i) + "," + String.valueOf(k+i));
                button.speciesIndex1 = i;
                button.speciesIndex2 = k+i;
                potButtons[k + potButtonsSpacer] = button;
                button.setBackground(Color.lightGray);
                button.addActionListener(buttonListener);
                gbl.setConstraints(button,gbc);
                leftPanePanel.add(button);
                gbc.gridx++;
            }
            potButtonsSpacer += getNumOSpecies()-i;
            gbc.gridx = 0;
            gbc.gridy++;
        }// end of button and label addition
           
        // As long as at least one species exists, the "Add," "Start," and "Remove" buttons are added
        if (getNumOSpecies() != 0) 
            addButtons();
        
        leftPanelHeight += 2*jButtonHeight;// Accounts for "Add," "Start," "Remove," and "Property Sheet" Buttons
        leftPanePanel.setMinimumSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
	    leftPanePanel.setPreferredSize(new java.awt.Dimension(leftPanelWidth, leftPanelHeight));
    }
}