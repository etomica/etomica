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
    public static final ActionListener OPEN = new OpenAction();
    
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
            
	        final SpaceSelectionFrame spaceFrame = new SpaceSelectionFrame(
                new SpaceSelectionFrame.Listener(){
                    public void spaceSelectionAction(Space s){
                        Simulation simulation = new Simulation(s);
                        Etomica.addSimulation(simulation);
                    }// end of spaceSelectionAction
                });// end of new Listener
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
    private static class OpenAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            FileDialog fd = new FileDialog(Etomica.DesktopFrame.etomicaFrame, "Open saved Simulation", FileDialog.LOAD);
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
    }// end of OpenAction
    
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
            SimulateActions.getSimulationEditorFrame().getSimulationEditor().resetAllComponentLists();
  //          Simulation.elementCoordinator.completed = false;
            Simulation.instance.speciesList = new LinkedList();
            Simulation.instance.potential1List = new LinkedList();
            Simulation.instance.potential2List = new LinkedList();
            Simulation.instance.integratorList = new LinkedList();
            Simulation.instance.phaseList = new LinkedList();
            Simulation.instance.controllerList = new LinkedList();
            Simulation.instance.displayList = new LinkedList();
            Simulation.instance.meterList = new LinkedList();
            Simulation.instance.deviceList = new LinkedList();
            Simulation.instance.graphicalElementList = new LinkedList();
//            Species.count = 0;
//            Simulation.instance = new simulate.Simulation(new simulate.Space2D());
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
    