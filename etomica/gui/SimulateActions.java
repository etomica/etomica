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
     * Static action listener for select the space of a simulation
     */
    public static final ActionListener SELECTSPACE = new SelectSpaceAction();
    
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
     * Static class that handles the select space action
     */
    public static class SelectSpaceAction implements ActionListener, java.io.Serializable {
        
        public void actionPerformed(ActionEvent event) {
	        final JInternalFrame spaceFrame = new JInternalFrame("Space",
	            true, // resizable
	            true, // closable
	            true, // maximizable
	            true); // iconifiable
    	        
            //BasicInternalFrameUI bifUI = new BasicInternalFrameUI(spaceFrame);
            //BasicToolBarUI btbUI = new BasicToolBarUI();
		    spaceFrame.getContentPane().setLayout(new GridLayout(0,1));
		    //((JComponent)spaceFrame).setUI(btbUI.createUI(((JComponent)spaceFrame)));
		    ButtonGroup dimension = new ButtonGroup();
		    
		    /**
		     * This section creates the radio buttons for all the classes that subclass space.class, makes
		     * them mutually exclusive, and adds buttonlisteners so that the selected button is known
		     */
		    for(int i=0; i<spaceClasses.length; i++) {
                String name = spaceClasses[i].getName();
                int idx = 13;//name.indexOf("r");  //strip off simulate.Space prefix
                name = name.substring(idx+1);
                JRadioButton button = new JRadioButton(name,i==0);
                dimension.add(button);
                spaceFrame.getContentPane().add(button);
            }// end of radio button creation
            
            /**
             * This section creates the "OK" button that when pressed creates an instance of 
             * simulation.instance, adds it to an internal frame, and opens the corresponding editor 
             * window. 
             */
            JButton oK = new JButton("OK");
            oK.addActionListener(new ActionListener(){
                    public void actionPerformed(ActionEvent e){
                        try{
                            EtomicaMenuBar.editSimulationItem.setEnabled(true);
                            spaceFrame.setClosed(true);
                            java.awt.Component[] allComponents = spaceFrame.getContentPane().getComponents();
                            // Run through all the space radio buttons
                            for(int j = 0; j < allComponents.length-1; j++){
                                // See which one is selected
                                if (((javax.swing.JRadioButton)allComponents[j]).isSelected()){
                                // Create SimulationFrame
                                    setApplet(new SimulationFrame(((javax.swing.AbstractButton)allComponents[j]).getText()));
                                    getApplet().addInternalFrameListener(new MyInternalFrameAdapter(){
                                        // When the simulation is closed, all of the tabs in the simulation editor
                                        // pane are reset to get ready for a new simulation
                                        public void internalFrameClosed(javax.swing.event.InternalFrameEvent evt){
                                            EtomicaMenuBar.selectSpaceItem.setEnabled(true);
                                            EtomicaMenuBar.editSimulationItem.setEnabled(false);
                                            SimEditorTabMenu.resetAllComponentLists();
                                            SimEditorTabMenu.setAllStart(false);
                                            SimEditorTabMenu.setAllRemove(false);
                                            SpeciesPotentialLinkPane.setNumOSpecies(0);
                                            SimEditorTabMenu.potential1Editor.update();
                                            SimEditorTabMenu.potential2Editor.update();
                                            try {
                                                getSimulationEditorFrame().setClosed(true);
                                            }
                                            catch(java.beans.PropertyVetoException pve){}
                                        }});// end of resetting simulation editor pane
                                    getApplet().reshape(520, 60, 470, 600);
                                    getApplet().setVisible(false);
                                    Etomica.DesktopFrame.desktop.add(getApplet());
                                // End creation of SimulationFrame    
                                    
                                // Create SimulationEditorFrame
                                    setSimulationEditorFrame(new SimulationEditorFrame());
                                    getSimulationEditorFrame().addInternalFrameListener(new MyInternalFrameAdapter(){
                                        public void internalFrameClosed(javax.swing.event.InternalFrameEvent evt){
                                            EtomicaMenuBar.editSimulationItem.setEnabled(true);
                                        }});
                                    getSimulationEditorFrame().setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
                                    getSimulationEditorFrame().reshape(10, 60, 500, 600);
                                    getSimulationEditorFrame().setVisible(true);
                                    Etomica.DesktopFrame.desktop.add(getSimulationEditorFrame());
                                // End creation of SimulationEditorFrame
                                
                                // Update MenuBars
                                    EtomicaMenuBar.editSimulationItem.setEnabled(false);
                                    EtomicaMenuBar.serEditItem.setEnabled(true);
                                    EtomicaMenuBar.serAppletItem.setEnabled(true);
                                    try{ 
                                        getSimulationEditorFrame().setSelected(true);
                                    }
                                    catch(PropertyVetoException exc){} // attempt was vetoed
                                }// end of manipulations necessary based on the selected space class
                            }// end of loop that runs through all of the space classes
                        }
                        catch(PropertyVetoException pve){}
                    }// end of actionPerformed
                    });// end of OK Button
            spaceFrame.getContentPane().add(oK);
		    spaceFrame.reshape(412, 200, 200, 200);
		    spaceFrame.setVisible(true);
		    Etomica.DesktopFrame.desktop.add(spaceFrame);

		    try{
		        spaceFrame.setSelected(true);
		    }
		    catch(PropertyVetoException e){} // attempt was vetoed
		}// end of actionPerformed
	
	    public static class MyInternalFrameAdapter extends javax.swing.event.InternalFrameAdapter implements java.io.Serializable {}
    }// end of SelectSpaceAction class

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
            catch(PropertyVetoException pve){}
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
    