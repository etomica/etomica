/**
 * FileActions
 *
 * The FileActions class creates static action listeners to the File drop-down menu of the 
 * EtomicaMenuBar.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */
package simulate.gui;

import simulate.*;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.beans.beancontext.BeanContextServicesSupport;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;
import javax.swing.JApplet;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDesktopPane;
import javax.swing.JFrame;
import javax.swing.ButtonGroup;
import javax.swing.JInternalFrame;
import javax.swing.JRadioButton;
import java.beans.PropertyVetoException;
import java.util.Hashtable;
import java.util.LinkedList;

public class FileActions {
    /**
     * Static action listener for creating a new simulation instance
     */
    public static final ActionListener NEW_SIMULATION = new NewSimulationAction();
    /**
     * Static action listener for serializing a simulation
     */
    public static final ActionListener SEREDIT = new SerEditAction();
    
    /**
     * Static action listener for serializing a simulation
     */
    public static final ActionListener SERAPPLET = new SerAppletAction();
    
    /**
     * Static action listener for loading a serialized component
     */
    public static final ActionListener LOAD = new LoadAction();
    
    /**
     * Static action listener for printing a simulation window
     */
    public static final ActionListener PRINT = new PrintAction();
    
    /**
     * Static action listener for clearing a simulation
     */
    public static final ActionListener CLEAR = new ClearAction();
    
    /**
     * Static action listener for exiting the etomica environment
     */
    public static final ActionListener EXIT = new ExitAction();
    
    /**
     * Static default file name for a serialized component
     */
    public final static String defaultSerializeAppletFile = "SimulationApplet.ser";

    /**
     * Static default file name for a serialized component
     */
    public final static String defaultSerializeEditFile = "EditableSimulation.ser";

    /**
     * Static default file name for the stored file
     */
    private final static String defaultStoreFile = "EditableSimulation.ser";
    
    transient static BeanContextServicesSupport bcss = new BeanContextServicesSupport();
    private static Hashtable serNames = new Hashtable();
    
    /**
     * Static class that creates a new simulation instance.
     */
    private static class NewSimulationAction implements ActionListener {
        
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
		    for(int i=0; i<SimulateActions.spaceClasses.length; i++) {
                String name = SimulateActions.spaceClasses[i].getName();
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
                                    SimulationFrame simulationFrame = new SimulationFrame(((javax.swing.AbstractButton)allComponents[j]).getText());
                                    Etomica.addSimulation(simulationFrame.simulation());
                                    SimulateActions.setApplet(simulationFrame);
                                    SimulateActions.getApplet().addInternalFrameListener(new SimulateActions.SelectSpaceAction.MyInternalFrameAdapter(){
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
                                                SimulateActions.getSimulationEditorFrame().setClosed(true);
                                            }
                                            catch(java.beans.PropertyVetoException pve){}
                                        }});// end of resetting simulation editor pane
                                    SimulateActions.getApplet().reshape(520, 60, 470, 600);
                                    SimulateActions.getApplet().setVisible(false);
                                    Etomica.DesktopFrame.desktop.add(SimulateActions.getApplet());
                                // End creation of SimulationFrame    
                                    
                                // Create SimulationEditorFrame
                                    SimulateActions.setSimulationEditorFrame(new SimulationEditorFrame());
                                    SimulateActions.getSimulationEditorFrame().addInternalFrameListener(new SimulateActions.SelectSpaceAction.MyInternalFrameAdapter(){
                                        public void internalFrameClosed(javax.swing.event.InternalFrameEvent evt){
                                            EtomicaMenuBar.editSimulationItem.setEnabled(true);
                                        }});
                                    SimulateActions.getSimulationEditorFrame().setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
                                    SimulateActions.getSimulationEditorFrame().reshape(10, 60, 500, 600);
                                    SimulateActions.getSimulationEditorFrame().setVisible(true);
                                    Etomica.DesktopFrame.desktop.add(SimulateActions.getSimulationEditorFrame());
                                // End creation of SimulationEditorFrame
                                
                                // Update MenuBars
                                    EtomicaMenuBar.editSimulationItem.setEnabled(false);
                                    EtomicaMenuBar.serEditItem.setEnabled(true);
                                    EtomicaMenuBar.serAppletItem.setEnabled(true);
                                    try{ 
                                        SimulateActions.getSimulationEditorFrame().setSelected(true);
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
        }//end of actionPerformed method
    }//end of NewSimulationAction class
    
    /**
     * Static class that handles the serialization action (Makes an editable version for the Etomica
     * environment
     */
    private static class SerEditAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
    	    FileDialog fd = new FileDialog(Etomica.DesktopFrame.etomicaFrame, "Serialize Component into File", FileDialog.SAVE);
	        // needed for a bug under Solaris...
	        fd.setDirectory(System.getProperty("user.dir"));
	        fd.setFile(defaultSerializeEditFile);
	        fd.show();
	        String fname = fd.getFile();
	        if (fname == null) {
	            return;
	        }
	        String dname = fd.getDirectory();
	        File file = new File(dname, fname);

	        try {
	            FileOutputStream f = new FileOutputStream(file);
	            ObjectOutputStream oos = new ObjectOutputStream(f);
	            // Ask the Wrapper to serialize the "naked" bean.
	            // as in copy()
	            oos.writeObject(SimulateActions.getApplet());//Simulation.instance);
	            oos.writeObject(SimulateActions.getSimulationEditorFrame());
	            oos.close();
	        } 
	        catch (Exception ex) {
	            error("Serialization of Component failed", ex);
	        }
        }// end of actionPerformed
        
        void error(String message, Throwable th) {
	        String mess = message + ":\n" + th;
	        System.err.println(message);
	        th.printStackTrace();

	        // Popup an ErrorDialog with the given error message.
	        new ErrorDialog(Etomica.DesktopFrame.etomicaFrame, mess);
        }// end of error
    }// end of SerEditAction
    
    /**
     * Static class that handles the serialization action (Makes an editable version for the Etomica
     * environment
     */
    private static class SerAppletAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            //MakeAppDlg appletDialog = new MakeAppDlg(Etomica.DesktopFrame.getTopBox(),((javax.swing.JPanel)((javax.swing.JInternalFrame)SimulateActions.getApplet()).getContentPane()));
    	    FileDialog fd = new FileDialog(Etomica.DesktopFrame.etomicaFrame, "Serialize Component into File", FileDialog.SAVE);
	        // needed for a bug under Solaris...
	        fd.setDirectory(System.getProperty("user.dir"));
	        fd.setFile(defaultSerializeAppletFile);
	        fd.show();
	        String fname = fd.getFile();
	        if (fname == null) {
	            return;
	        }
	        String dname = fd.getDirectory();
	        File file = new File(dname, fname);

	        try {
	            FileOutputStream f = new FileOutputStream(file);
	            ObjectOutputStream oos = new ObjectOutputStream(f);
	            // Ask the Wrapper to serialize the "naked" bean.
	            // as in copy()
	            JApplet applet = new JApplet();
	            applet.getContentPane().add(Simulation.instance);
	            oos.writeObject(applet);//Simulation.instance);
	            oos.close();
	        } 
	        catch (Exception ex) {
	            error("Serialization of Component failed", ex);
	        }
        }// end of actionPerformed
        
        void error(String message, Throwable th) {
	        String mess = message + ":\n" + th;
	        System.err.println(message);
	        th.printStackTrace();

	        // Popup an ErrorDialog with the given error message.
	        new ErrorDialog(Etomica.DesktopFrame.etomicaFrame, mess);
        }// end of error
    }// end of SerAppletAction
    
    /**
     * Static class that handles the load a serialized component action
     */
    private static class LoadAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            FileDialog fd = new FileDialog(Etomica.DesktopFrame.etomicaFrame, "Load saved Simulation", FileDialog.LOAD);
	        // needed for a bug under Solaris...
	        fd.setDirectory(System.getProperty("user.dir"));
	        fd.setFile(defaultStoreFile);
	        fd.show();
	        String fname = fd.getFile();
	        if (fname == null) {
	            return;
	        }
	        String dname = fd.getDirectory();
	        File file = new File(dname, fname);

	        // OK, loaded all the classes -- now, instantiate them...
	        try {
	            // get the one object as a bean...
	            FileInputStream f = new FileInputStream(file);
	            ObjectInputStream ois = new ObjectInputStream(f);
                Object applet = ois.readObject();
                Object editor = ois.readObject();
                ois.close();
               
                SimulateActions.setApplet((SimulationFrame)applet);
                Etomica.DesktopFrame.desktop.add(SimulateActions.getApplet());
                Object instance = ((JInternalFrame)SimulateActions.getApplet()).getContentPane().getComponent(0);

                SimulateActions.setSimulationEditorFrame((SimulationEditorFrame)editor);
                Etomica.DesktopFrame.desktop.add(SimulateActions.getSimulationEditorFrame());
                EtomicaMenuBar.editSimulationItem.setEnabled(false);
                EtomicaMenuBar.serEditItem.setEnabled(true);
                EtomicaMenuBar.serAppletItem.setEnabled(true);
	        } 
	        catch(java.io.IOException ioe){ ioe.printStackTrace(); }
	        catch(ClassNotFoundException cnfe){ cnfe.printStackTrace(); }
        }// end of actionPerformed
    }// end of LoadAction
    
    /**
     * Static class that handles the print a component action
     */
    private static class PrintAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            java.awt.PrintJob pj = java.awt.Toolkit.getDefaultToolkit().getPrintJob(Etomica.DesktopFrame.etomicaFrame, "Printing Test", (java.util.Properties)null);

    	    if (pj != null) {
//	    	    Graphics g;
//		        Dimension pageDim = pj.getPageDimension();
//		        int pageRes = pj.getPageResolution();
//		        boolean lastFirst = pj.lastPageFirst();

		        Graphics gr  = pj.getGraphics();
    		
		        if (gr!=null) {
		            // We print all components one after the other in the
		            // same page. 		  
		            int count = Etomica.DesktopFrame.etomicaFrame.getComponentCount();
		            System.out.println(count);
		            for (int i=0; i<count; i++) {
		                Wrapper wr = new Wrapper(SimulateActions.getApplet(),"SimulationFrame1","simulate.gui.SimulationFrame");
		                Object bean = wr.getBean();
		                if (bean instanceof Component) {
		                    Component c = (Component) java.beans.Beans.getInstanceOf(bean, Component.class);
		                    Dimension d = c.getSize();
		                    java.awt.Point o = wr.getLocation();
		                    java.awt.Image offScreen = c.createImage(d.width, d.height);
		                    c.paint(offScreen.getGraphics());
		                    gr.drawImage(offScreen, o.x, o.y, java.awt.Color.white, null);
		                }
		            }
		        }
		        else System.err.println("Could not get Graphics handle.");
        			     
		        gr.dispose();
	            pj.end();
	        } 
	        else System.err.println("PrintJob cancelled.");     

        }// end of actionPerformed
    }// end of PrintAction
    
    /**
     * Static class that handles the clear the simulation action
     */
    private static class ClearAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            try {
                SimulateActions.getApplet().setClosed(true);
                SimulateActions.getSimulationEditorFrame().setClosed(true);
            }
            catch (java.beans.PropertyVetoException pve){}
            SimulateActions.getApplet().getContentPane().removeAll();
            SimEditorTabMenu.resetAllComponentLists();
  //          Simulation.elementCoordinator.completed = false;
            Simulation.speciesList = new LinkedList();
            Simulation.potential1List = new LinkedList();
            Simulation.potential2List = new LinkedList();
            Simulation.integratorList = new LinkedList();
            Simulation.phaseList = new LinkedList();
            Simulation.controllerList = new LinkedList();
            Simulation.displayList = new LinkedList();
            Simulation.meterList = new LinkedList();
            Simulation.deviceList = new LinkedList();
            Simulation.graphicalElementList = new LinkedList();
//            Species.count = 0;
            Simulation.instance = new simulate.Simulation(new simulate.Space2D());
            //SimulateActions.getApplet().getContentPane().add(Simulation.instance);
        }// end of actionPerformed
    }// end of ClearAction
    
    /**
     * Static class that handles the exit the Etomica environment action
     */
    private static class ExitAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            System.exit(0);
        }// end of actionPerformed
    }// end of ExitAction
}// end of FileActions class
    