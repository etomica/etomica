/**
 * DefineAtomPotentialFrame
 *
 * The DefineAtomPotentialFrame class is a splitpane that lists all the simulation components of a 
 * respective potential category (potential 1 or potential 2) on the leftside, and all of the added 
 * components from the leftside list on the rightside in a JList format.
 *
 * @author Bryan C. Mihalick
 * 9/21/00
 */

package simulate.gui;

import simulate.*;
import java.awt.Color;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JRadioButton;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import java.io.File;
import java.io.FilenameFilter;

public class DefinePotentialFrame extends JInternalFrame {
    /**
     * Lists all of the simulation components corresponding to the respective tabs name.  These are listed
     * as radio buttons.
     */
    final javax.swing.JPanel atomButtonPanel = new javax.swing.JPanel();
    
    /**
     * Makes it possible to determine which radio button was selected when the "add" button is pressed
     */
    final ButtonListener buttonListener = new ButtonListener();
    
    /**
     * Holds all constraints needed for displaying the next awt or swing component
     */
    final GridBagConstraints gbc = new GridBagConstraints();
    
    /**
     * Determines how to display an awt or swing component by using gbc from above
     */
    final GridBagLayout gbl = new GridBagLayout();

    /**
     * Button for starting the simulation.  Only enabled if a sufficient number of simulation components
     * have been added to make a working simulation
     */
    final JButton done = new JButton("Done");
    
    /**
     * Index of the potential class that the currently selected JButton in the left pane corresponds to
     */
    int currentIndex = 0;
    
    /**
     * Index of object that is currently selected in the JList on the right pane
     */
    int currentSelection = -1;
    
    /**
     * The title that will be displayed on the property sheet when a component is selected
     */
    String title;
    
    /**
     * If true, the current simulation component has already been added
     */
    boolean added = false;
    
    /**
     * Frame containing radioButtons for all of the Potential classes
     */
    public PotentialFrame potentialFrame;
    
    /**
     * If true, at least one species was added to the simulation.  This makes sure that the
     * actionListener of the remove button is only added once.
     */
    boolean firstSpecies = true;
    
    /**
     * Determines if Potential Arrays need to be instantiated or not
     */
    public static boolean makePotArrays = true;
    
    /**
     * Static array of all simulation components that extend potential1.class
     */
    public static Class[] potentialClasses;
    
    /**
     * Handle to the species objects that correpsonds to the currently selected JButton on the 
     * SpeciesPotentialLinkPane
     */
    public static Species species1, species2;
    
    /**
     * Indices of the species objects that correpsonds to the currently selected JButton on the 
     * SpeciesPotentialLinkPane
     */
    public static int speciesIndex1, speciesIndex2;
    
    public static JButton currentButton;
    JButton remove = new JButton();
    private static int P1IDnumber = 0, P2IDnumber = 0;
    
    protected SimulationEditor simulationEditor;

    /**
     * Constructor that creates all of the left pane's radiobuttons and JButtons, as well as, the right 
     * pane's scrollpane and JList.  It also creates all the listeners for these swing components so that
     * the simulation can be updated as needed.
     */
    public DefinePotentialFrame(SimulationEditor ed){
        simulationEditor = ed;
        setBounds(515,490,PotentialFrame.atomPairPotArray[0].length*50 + 150,PotentialFrame.atomPairPotArray.length*25 + 100);
        setResizable(true);
        setVisible(true);
        gbc.gridwidth = gbc.REMAINDER;
        gbc.gridx = 0;
        gbc.gridy = 0;
        
        /**
         * This section sets up the JList of the right pane.  When an object from the list is selected,
         * its properties are displayed in the property sheet. 
         */
        getContentPane().add(atomButtonPanel);
        createGrid();
    }

    /**
     * This method is called when a species is added or removed from the species pane.  It figures out
     * the number of possible potential interactions, creates buttons for each of these interactions,
     * determines the necessary layout of the buttons, and displays them.
     */
    public void createGrid(){
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 1;
        gbc.anchor = gbc.WEST;
        atomButtonPanel.removeAll();
        atomButtonPanel.setLayout(gbl);
        JLabel label1 = new JLabel("   ");
        atomButtonPanel.add(label1);
        gbc.gridx++;
        
        int potButtonsSpacer = 0;
        
        for (int i = 0; i < species2.getAtomsPerMolecule(); i++){
            JLabel label2 = new JLabel("Atom " + String.valueOf(i));
            gbl.setConstraints(label2,gbc);
            atomButtonPanel.add(label2);
            gbc.gridx++;
        }

        gbc.gridx = 0;
        gbc.gridy++;
        for (int i = 0; i < species1.getAtomsPerMolecule(); i++){
            JLabel label3 = new JLabel("Atom " + String.valueOf(i));
            gbl.setConstraints(label3,gbc);
            atomButtonPanel.add(label3);
            gbc.gridx++;
            // add buttons with listeners and change color
            for (int k = 0; k < species2.getAtomsPerMolecule(); k++){
                SpeciesPairButton button = new SpeciesPairButton(String.valueOf(i) + "," + String.valueOf(k));
                button.speciesIndex1 = i;
                button.speciesIndex2 = k;
                button.setBackground(Color.lightGray);
                button.addActionListener(buttonListener);
                gbl.setConstraints(button,gbc);
                atomButtonPanel.add(button);
                gbc.gridx++;
            }
            gbc.gridx = 0;
            gbc.gridy++;
        }// end of button and label addition
        
        /**
         * This section adds the start button, but it's disabled until enough simulation components
         * have been added to make a working simulation.
         */
        gbc.gridx = 0;
        gbc.gridy++;
	    done.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e){
	            Potential2 p2IdealGas = new P2IdealGas();
	            for (int i = 0; i < PotentialFrame.atomPairPotArray.length; i++){
	                for (int j = 0; j < PotentialFrame.atomPairPotArray[0].length; j++){
	                    if (PotentialFrame.atomPairPotArray[i][j] == null)
	                        PotentialFrame.atomPairPotArray[i][j] = (Class)p2IdealGas.getClass();
	                }
	            }
	            if (PotentialFrame.getPotentialEditor().getTitle() == "Potential1") {
	                P1DefinedPotential newPotential = new P1DefinedPotential(PotentialFrame.atomPairPotArray);
	                newPotential.setName("Bryan's P1 Potential - " + Integer.toString(P1IDnumber++));
	                simulationEditor.potential1Editor.componentList.addElement(newPotential);
                    simulationEditor.getSimulation().potential1[speciesIndex1] = newPotential;
                }
                else {
	                P2DefinedPotential newPotential = new P2DefinedPotential(PotentialFrame.atomPairPotArray);
	                newPotential.setName("Bryan's P2 Potential - " + Integer.toString(P2IDnumber++));
                    simulationEditor.potential2Editor.componentList.addElement(newPotential);
                    simulationEditor.getSimulation().potential2[speciesIndex1][speciesIndex2] = newPotential;
                    simulationEditor.getSimulation().potential2[speciesIndex2][speciesIndex1] = newPotential;
                }
	            
	            
	            
	            try {
                    ((JInternalFrame)Etomica.DesktopFrame.desktop.getComponent(0)).setClosed(true);
                }
                catch(java.beans.PropertyVetoException pve){}            
            }});
        gbl.setConstraints(done,gbc);
        atomButtonPanel.add(done);
        // end of start button
    }

    /**
     * Sets index of first species that contains the atoms to have new potential applied upon
     */
    public static void setSpeciesIndex1(int s1){ speciesIndex1 = s1; }

    /**
     * Sets index of second species that contains the atoms to have new potential applied upon
     */
    public static void setSpeciesIndex2(int s2){ speciesIndex2 = s2; }
    
    /**
     * Sets first species that contains the atoms to have new potential applied upon
     */
    public static void setSpecies1(Species s){ species1 = s; }

    /**
     * Sets second species that contains the atoms to have new potential applied upon
     */
    public static void setSpecies2(Species s){ species2 = s; }
    
    /**
     * This class allows buttons to hold a handle to the potential class that is associated to the button,
     * as well as the corresponding index of that potential in the JList
     */
    public class SpeciesPairButton extends javax.swing.JButton{
        Class potential = null;
        int index = -1;
        int speciesIndex1 = -1;
        int speciesIndex2 = -1;
        
        SpeciesPairButton(String name){
            super(name);
        }// end of SpeciesPairButton constructor        
    }// end of SpeciesPairButton class
    
    /**
     * Extension class of ActionListener that saves a handle to the buttons that are currently 
     * selected.  It also turns the selected buttons to red, and then instantiates a PotentialFrame
     * window that lists all of the possible potential classes that can be added to the simulation.
     */
	private class ButtonListener implements ActionListener {
	    
	    public void actionPerformed(ActionEvent evt) {
	        currentButton = ((JButton)evt.getSource());
	        ((SpeciesPairButton)currentButton).index = currentIndex;
	        currentSelection = -1;
	        ((java.awt.Component)evt.getSource()).setBackground(Color.red);

            if (makePotArrays) {
//                simulationEditor.getSimulation().potential1 = new Potential1[SpeciesPotentialLinkPane.getNumOSpecies()];
//                simulationEditor.getSimulation().potential2 = new Potential2[SpeciesPotentialLinkPane.getNumOSpecies()][SpeciesPotentialLinkPane.getNumOSpecies()];
                makePotArrays = false;
            }
	        potentialFrame = new PotentialFrame(simulationEditor);
	        potentialFrame.setTitle(getTitle());
            Etomica.DesktopFrame.desktop.add(potentialFrame);
	        try {
	            potentialFrame.setSelected(true);
	        }
	        catch(java.beans.PropertyVetoException pve){}
        }// end of actionPerformed
	}// end of ButtonListener class
}// end of DefineAtomPotentialFrame class