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
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import java.io.File;
import java.io.FilenameFilter;

public abstract class PotentialFrame extends javax.swing.JInternalFrame {
    /**
     * JPanel used to hold all of the components
     */
    JPanel potentialPanel = new JPanel();
    
    /**
     * Scrollpane that holds the potentialPanel;
     */
    JScrollPane scroller = new JScrollPane(potentialPanel);
    
    /**
     * Button for adding objects to the JList of the right pane.  
     */
    final JButton addToSim = new JButton("Add");
    
    /**
     * Button for canceling out of the PotentialFrame
     */
    final JButton cancel = new JButton("Cancel");
    
    /**
     * Handle to the 'Remove' button from the SpeciesEditorLinkPane
     */
    public static JButton removeButton = new JButton();
    
    /**
     * Set to true if added potential should be instantiated.  Will be false if added potential is part
     * of a serialized simulation, and will be true otherwise.
     */
    boolean instantiate = true;

    /**
     * Array of class names created in the SimulateActions classe's static block
     */
    Class[] className;
    
    /**
     * Array of Atom-Atom Potential Interactions
     */
    public static Class[][] atomPairPotArray; 
        
    /**
     * ButtonGroup that allows only one of the JRadioButtons to be selected at a time
     */
    final ButtonGroup simComponents = new ButtonGroup();
        
    /**
     * ButtonListener that keeps track of the currently selected MyRadioButton
     */
    final ButtonListener buttonListener = new ButtonListener();
        
    /**
     * Holds that current constraints used determining the position of the added component
     */
    final GridBagConstraints gbc = new GridBagConstraints();
        
    /**
     * Handle to the current GridBagLayout
     */
    final GridBagLayout gbl = new GridBagLayout();
        
    /**
     * Handle to the MyRadioButton that is currently selected
     */
    MyRadioButton currentButton;
        
    /**
     * Handle to the most recently instantiated potential class
     */
    Object component;
    
    /**
     * Handle to the corresponding PotentialEditorPane
     */
    protected static SpeciesPotentialLinkPane potentialEditor;
    
    /**
     * Height of the potentialPanel above
     */
    protected int panelHeight = 0;
    
    /**
     * Width of the potentialPanel above
     */
    protected int panelWidth = 175;
    
    /**
     * Height of one JRadioButton.  Used in determining the panelHeight of the potentialPanel.
     */
    static final int radioButtonHeight = 25;
    
    /**
     * Height of one JButton.  Used in determining the panelHeight of the potentialPanel.
     */
    static final int jButtonHeight = 35;
    
    /**
     * Static array of all simulation components that extend potential.class
     */
    public static Class[] potentialClasses; 

    /**
     * Static array of all simulation components that extend potential.class
     */
    public static Class[] potential1Classes; 
    
    /**
     * Static array of all simulation components that extend potential.class
     */
    public static Class[] potential2Classes; 

    /**
     * Constructor that takes the potential type (P1 or P2) and the removeButton from the Potential tab
     * of the SimEditorTabMenu that called this constructor.  It determines the necessary size based on
     * the type of potentials being used, sets the title, and adds the JInternalFrame to the JDesktopPane
     */
    public PotentialFrame(){
        setTitle("Potentials");  
        setResizable(true);
        
        gbc.gridwidth = 3;
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.anchor = gbc.WEST;

		potentialPanel.setLayout(gbl);
        getContentPane().add(scroller);
    }// end of PotentialFrame constructor

    public static void setRemoveButton(JButton b) { removeButton = b; }
    
    public void makeRadioButtons(Class[] className, int idx) {
		/**
		 * This section creates all of the left pane's radio buttons, adds listeners to these buttons,
		 * and adds them to the pane in a grid bag layout.
		 */
		for(int i=0; i < className.length; i++) {
            String name = className[i].getName();
            name = name.substring(idx+1);            
            MyRadioButton button = new MyRadioButton(name,false,className[i]);
            button.addActionListener(buttonListener);
            gbl.setConstraints(button, gbc);
            simComponents.add(button);
            potentialPanel.add(button);
            gbc.gridy++;
            panelHeight += radioButtonHeight;// Accounts for height of the new radioButton
        }// end radio button addition
    }

    /**
     * This section creates the "Add" button which when pressed creates an instance of the class that
     * corresponds to the currently selected radio button of the left pane.  It then adds this object
     * to the JList of the rightpane and the simulation.instance object.
     */ 
    public final void addAddButton(){
        //gbc.gridx = 0;
        gbc.gridwidth = 1;
        gbl.setConstraints(addToSim, gbc);
        potentialPanel.add(addToSim);
    }// end of add button

    public final void addCancelButton(){
        gbc.gridx++;
        gbl.setConstraints(cancel, gbc);
        potentialPanel.add(cancel);
    }// end of cancel button

    public static void setPotentialEditor(SpeciesPotentialLinkPane p) { potentialEditor = p; }
    public static SpeciesPotentialLinkPane getPotentialEditor() { return potentialEditor; }
    
    /**
     * Class that has all the functionality of a JRadioButton with the addition of being able to
     * remember the class that is associated with that JRadioButton
     */
    protected class MyRadioButton extends JRadioButton {
        Object cls;
            
        MyRadioButton(String label, boolean selected, Object cls) {
            super(label,selected);
            this.cls = cls;
        }
    }// end of MyRadioButton
        
    /**
     * Class that allows each MyRadioButton to update the currently selected button when they are selected
     */
	protected class ButtonListener implements ActionListener, java.io.Serializable {
    	    
	    public void actionPerformed(ActionEvent evt) {
	        currentButton = ((MyRadioButton)evt.getSource());
        }//end of actionPerformed
	}//end of ButtonListener
	
	public class MyActionListener implements java.awt.event.ActionListener, java.io.Serializable {
	    public void actionPerformed(java.awt.event.ActionEvent ae){}
	}
	
	static {
    	// Initialization of potentialClasses array
	    File dir = new File(simulate.Default.CLASS_DIRECTORY);
	    String[] files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Potential")
                    && name.endsWith("class")
                    && !name.startsWith("PotentialField$")
                    && !name.startsWith("PotentialFieldGravity$")
                    && !name.startsWith("PotentialRoughSphere$")
                    && !name.startsWith("Potential1")
                    && !name.startsWith("Potential2")
	                && !name.startsWith("Potential$")
                    && !name.startsWith("Potential.class");
            }});
        potentialClasses = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        potentialClasses[i] = null;
	        try{
	            potentialClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of potentialClasses array

    	// Initialization of potential1Classes array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("P1")
                    && name.endsWith("class")
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && name.indexOf("$") == -1;}
	        });
        potential1Classes = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        potential1Classes[i] = null;
	        try{
	            potential1Classes[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of potential1Classes array
    	
    	// Initialization of potential2Classes array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Potential")
	                && name.endsWith("class")
                    && !name.startsWith("PotentialField")
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
                    && !name.startsWith("Potential1")
                    && !name.startsWith("Potential2")
	                && !name.startsWith("Potential$")
                    && !name.startsWith("Potential.class")
	                && name.indexOf("$") == -1;}
	        });
        potential2Classes = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        potential2Classes[i] = null;
	        try{
	            potential2Classes[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of potential2Classes array
   	}// end of static block used for potential class array initialization
}// end of PotentialFrame class