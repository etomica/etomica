/**
 * PreferenceFrame.java
 *
 * The PreferenceFrame is a JInternalFrame that lists the current setup preferences for the
 * Etomica environment.  It contains a JTabbedPane that has a different pane for each property
 * type. 
 *
 * @author Bryan C. Mihalick
 *
 * 1/17/01
 */
 
package etomica.gui;

import etomica.*;
import etomica.units.UnitSystem;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

public class PreferenceFrame extends javax.swing.JInternalFrame {
    public static JTabbedPane tabbedPane = new JTabbedPane();
    public static JPanel unitPanel = new JPanel();
    public static Class[] unitSystemClasses;
    private static final ButtonGroup buttonGroup = new ButtonGroup();
    private static final JButton ok = new JButton("ok");
    private static ButtonListener buttonListener = null;
    private static MyButton currentButton = null;
    /**
     * Holds all constraints needed for displaying the next awt or swing component
     */
    final GridBagConstraints gbc = new GridBagConstraints();
    
    /**
     * Determines how to display an awt or swing component by using gbc from above
     */
    final GridBagLayout gbl = new GridBagLayout();

    PreferenceFrame(){
        super("Preferences",true,false,false,true);//resizable and iconifiable
        setBounds(412,288,200,200);
        unitPanel.setSize(200, 200);
        getContentPane().setLayout(gbl);
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 3;
        gbc.anchor = gbc.WEST;
        buttonListener = new ButtonListener();
        
	    // Creates a radiobutton for each type of UnitSystem
	    for(int i=0; i<unitSystemClasses.length; i++) {
            String name = unitSystemClasses[i].getName();
            int idx = 24;//strip off etomica.units.UnitSystem$ prefix
            name = name.substring(idx+1);
            MyButton button = new MyButton(unitSystemClasses[i],name,i==0);
            button.addActionListener(buttonListener);
            gbl.setConstraints(button, gbc);
            buttonGroup.add(button);
            unitPanel.add(button);
            gbc.gridy++;
        }// end of radio button creation

        tabbedPane.addTab("UnitSystem", unitPanel);
        getContentPane().add(tabbedPane);
        
        ok.addActionListener(new java.awt.event.ActionListener(){
            public void actionPerformed(java.awt.event.ActionEvent e){
                try {
                    UnitSystem system = ((UnitSystem)((Class)currentButton.unitSystemClass).newInstance());
                    Simulation.instance.setUnitSystem(system);
                    setClosed(true);
                }
	            catch(InstantiationException exc) { exc.printStackTrace(); }
	            catch(IllegalAccessException exc) { exc.printStackTrace(); }
	            catch(java.beans.PropertyVetoException pve) {}
            }
        });
        gbc.gridx = 0;
        gbc.gridwidth = 1;
        gbl.setConstraints(ok, gbc);
        getContentPane().add(ok);
        
        setVisible(true);
        Etomica.DesktopFrame.desktop.add(this);
        try { setSelected(true);}
        catch (java.beans.PropertyVetoException pve){}
    }//end of PreferenceFrame constructor
    
    /**
     * Simple extension of JRadioButton to hold its corresponding UnitSystem class object.
     */
    private static class MyButton extends javax.swing.JRadioButton {
        Class unitSystemClass;
        MyButton(Class c, String name, boolean selected) {
            super(name, selected);
            unitSystemClass = c;
        }
    }

    /**
     * Extension class of ActionListener that saves a handle to the radio button that is currently selected
     */
	private class ButtonListener implements java.awt.event.ActionListener {
	    
	    public void actionPerformed(java.awt.event.ActionEvent evt) {
	        currentButton = ((MyButton)evt.getSource());
        }//end of actionPerformed
	}//end of ButtonListener

	static {
	    java.io.File dir = new java.io.File(etomica.Default.CLASS_DIRECTORY + "/units");
	    String[] files = dir.list(new java.io.FilenameFilter() {
	        public boolean accept(java.io.File d, String name) {
	                return name.startsWith("UnitSystem$")
	                && name.endsWith("class");}
	        });
	    unitSystemClasses = new Class[files.length];
	    for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        unitSystemClasses[i] = null;
	        try{
	            unitSystemClasses[i] = Class.forName("etomica.units."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End of initialization of UnitSystems array
	}// end of static block
}//end of PreferenceFrame class
