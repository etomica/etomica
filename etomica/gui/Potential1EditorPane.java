package simulate.gui;

import java.awt.Color;
import javax.swing.JButton;
import javax.swing.JLabel;

public class Potential1EditorPane extends SpeciesPotentialLinkPane {
    
    Potential1EditorPane(){
        super();
        setTitle("Potential1");
    }
    
    public void update(){
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 1;
        gbc.anchor = gbc.WEST;
        leftPanelHeight = (getNumOSpecies()+3)*jButtonHeight;
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
        currentButtons = new JButton[getNumOSpecies()+1];
        potButtons = new JButton[getNumOSpecies()+1];
 
        gbc.gridx = 0;
        gbc.gridy++;
        // add buttons and labels
        for (int i = 0; i < getNumOSpecies(); i++){
            JLabel label2 = new JLabel("Species " + String.valueOf(i));
            gbl.setConstraints(label2,gbc);
            leftPanePanel.add(label2);
            gbc.gridx++;

            SpeciesPairButton button = new SpeciesPairButton(String.valueOf(i));
            button.speciesIndex1 = i;
            potButtons[i] = button;
            button.setBackground(Color.lightGray);
            button.addActionListener(new ButtonListener());
            gbl.setConstraints(button,gbc);
            leftPanePanel.add(button);
            gbc.gridx = 0;
            gbc.gridy++;
        }// end of button and label addition
        
        // As long as at least one species exists, the "Add," "Start," and "Remove" buttons are added
        if (getNumOSpecies() != 0)
            addButtons();
    }// end of update method
}
