/**
 * SimulateActions
 *
 * The SimulateActions class creates static action listeners to the simulate drop-down menu of the
 * EtomicaMenuBar.
 *
 * @author Bryan C. Mihalick
 * 8/16/00
 */

package simulate.gui;

import java.awt.Component;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDesktopPane;
import javax.swing.JFrame;
import javax.swing.JRadioButton;
import javax.swing.JToolBar;
import javax.swing.JInternalFrame;
import java.beans.PropertyVetoException;
import java.io.File;
import java.io.FilenameFilter;
import simulate.*;

public class SimulateActions {
        
    /**
     * Static action listener for opening the editor window corresponding to the current simulation
     */
    public static final ActionListener EDITSIMULATION = new EditSimulationAction();
    
    /**
     * Static array of all simulation components that extend space.class
     */
    public static Class[] spaceClasses;
    
    /**
     * Static array of all simulation components that extend potential1.class
     */
    public static Class[] potential1Classes;
    
    /**
     * Static array of all simulation components that extend potential2.class
     */
    public static Class[] potential2Classes;
    
    /**
     * Static array of all simulation components that extend species.class
     */
    public static Class[] speciesClasses;
    
    /**
     * Static array of all simulation components that extend integrator.class
     */
    public static Class[] integratorClasses;
    
    /**
     * Static array of all simulation components that extend phase.class
     */
    public static Class[] phaseClasses;
    
    /**
     * Static array of all simulation components that extend controller.class
     */
    public static Class[] controllerClasses;
    
    /**
     * Static array of all simulation components that extend display.class
     */
    public static Class[] displayClasses;
    
    /**
     * Static array of all simulation components that extend meter.class
     */
    public static Class[] meterClasses;
    
    /**
     * Static array of all simulation components that extend device.class
     */
    public static Class[] deviceClasses;
    
    /**
     * Static handle to the internal frame that contains the simulation.instance object
     */
    private static JInternalFrame appletFrame;
    
    /**
     * Static handle to the editor window corresponding to the current simulation
     */
    private static SimulationEditorFrame editorFrame;

    /**
     * Static mutator method for setting the JInternalFrame that contains the simulation.instance object
     */
    public static void setApplet(SimulationFrame applet) {appletFrame = applet;}
    
    /**
     * Static accessor method that returns the JInternalFrame that contains the simulation.instance object
     */
    public static JInternalFrame getApplet(){ return appletFrame; }

    /**
     * Static mutator method that sets the editor window that corresponds to the current simulation
     */
    public static void setSimulationEditorFrame(SimulationEditorFrame frame){ editorFrame = frame; }
    
    /**
     * Static accessor method that gets the editor window that corresponds to the current simulation
     */
    public static SimulationEditorFrame getSimulationEditorFrame(){ return editorFrame; }

    /**
     * Static class that handles the edit simulation action
     */
    private static class EditSimulationAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            EtomicaMenuBar.editSimulationItem.setEnabled(false);
            try {
                SimulateActions.getSimulationEditorFrame().setClosed(false);
                SimulateActions.getSimulationEditorFrame().setSelected(true);
            }
            catch(java.beans.PropertyVetoException pve){}
            Etomica.DesktopFrame.desktop.add(SimulateActions.getSimulationEditorFrame());
        }// end of actionPerformed
    }// end of EditSimulationAction class
    
    /**
     * Class describing actions to be performed if a different simulation is selected.
     * Each instantiated simulation has its own instance of this action, which is called
     * if it becomes the current simulation.
     */
    public static class SelectSimulationAction implements ActionListener {
        Simulation simulation;
        /**
         * Constructor takes as argument the simulation associated with this action.
         */
        public SelectSimulationAction(Simulation sim) {
            simulation = sim;
        }
        public void actionPerformed(ActionEvent event) {
            Simulation.instance = simulation;
        }
    }

    /**
     * This static block sets up all of the static arrays that hold all of the simulation component
     * classes
     */
    static {
        // Initialization of spaceClasses array
        
        //spaceclasses is moved to Etomica, so if everything works this segment can be deleted
	    File dir = new File(simulate.Default.CLASS_DIRECTORY);
	    String[] files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
	                return name.startsWith("Space")
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && !name.equals("Space.class");}
	        });
	    spaceClasses = new Class[files.length];
	    for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        spaceClasses[i] = null;
	        try{
	            spaceClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End of initialization of spaceClasses array
    	
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
                return name.startsWith("P2")
	                && name.endsWith("class")
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
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
    	
    	// Initialization of speciesClasses array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Species")
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && !name.equals("Species.class");}
	        });
        speciesClasses = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        speciesClasses[i] = null;
	        try{
	            speciesClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of speciesClasses array
    	    
	    // Initialization of integratorClasses array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Integrator")
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && !name.equals("Integrator.class");}
	        });
        integratorClasses = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        integratorClasses[i] = null;
	        try{
	            integratorClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of integratorClasses array
    	    
	    // Initialization of phaseClasses array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Phase")
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && !name.startsWith("PhaseAction")
	                && !name.startsWith("PhaseEvent");}	                	                
	        });
        phaseClasses = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        phaseClasses[i] = null;
	        try{
	            phaseClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of phaseClasses array
    	    
	    // Initialization of controllerClasses array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Controller")
	                && name.endsWith("class")
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && name.indexOf("$") == -1;}
	        });
        controllerClasses = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        controllerClasses[i] = null;
	        try{
	            controllerClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of controllerClasses array
	    
	    // Initialization of displayClasses array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Display")
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && !name.startsWith("DisplayPhaseEvent")
	                && !name.startsWith("DisplayPhaseListener")
	                && !name.equals("Display.class");}
	        });
        displayClasses = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        displayClasses[i] = null;
	        try{
	            displayClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of displayClasses array
	    
	    // Initialization of meterClasses array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Meter")
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && !name.startsWith("MeterAbstract")
	                && !name.equals("Meter.class");}
	        });
        meterClasses = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        meterClasses[i] = null;
	        try{
	            meterClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of meterClasses array
	    
	    // Initialization of deviceClasses array
	    files = dir.list(new FilenameFilter() {
	        public boolean accept(File d, String name) {
                return name.startsWith("Device")
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class")
	                && !name.equals("Device.class");}
	        });
        deviceClasses = new Class[files.length];
        for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        deviceClasses[i] = null;
	        try{
	            deviceClasses[i] = Class.forName("simulate."+files[i]);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	    }// End initialization of deviceClasses array
	}// end of static block for class array initiliazation    
}// end of SimulateActions class
    