/**
 * ViewActions
 *
 * The ViewActions class is responsible for creating static action listeners to the View drop-down menu
 * of the EtomicaMenuBar.
 *
 * @author Bryan C. Mihalick
 * 8/15/00
 */

package etomica.gui;

import java.awt.GridLayout;
import java.awt.Frame;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JInternalFrame;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import java.beans.PropertyVetoException;
import etomica.*;

public class ViewActions {
    /**
     * Static action listener for the workbook event which puts all of the source files on a tabbedPane
     */
    public static final ActionListener WORKBOOK = new WorkbookAction();
    
    /**
     * Static action listener for the property list event which opens the property list 
     */
    public static final ActionListener PROPERTYLIST = new PropertyListAction();
    
    /**
     * Static action listener for the component library event which opens the component library
     */
    public static final ActionListener COMPONENTLIBRARY = new ComponentLibraryAction();
    
    /**
     * Static action listener for the class browser event which opens the class browser
     */
    public static final ActionListener CLASSBROWSER = new ClassBrowserAction();
    
    /**
     * Static action listener for the messages event which opens the messages window
     */
    public static final ActionListener MESSAGES = new MessagesAction();
    
    /**
     * Static action listener for events that open the simulation window.
     */
    public static final ActionListener SIMULATION = new SimulationAction();
        
    /**
     * Static class that handles the workbook event
     */
    private static class WorkbookAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of WorkbookAction class

    /**
     * Static class that handles the property list event
     */
    private static class PropertyListAction implements ActionListener {
        /**
         * handle to the property sheet of the currently selected object
         */
        public PropertySheet propertySheet; 
        boolean added = false;
        
        PropertyListAction(){
            propertySheet = new PropertySheet(null, 515, 60);
            Etomica.setPropertySheet(propertySheet);
            SimulationEditorPane.setPropertySheet(propertySheet);
            PotentialEditorPane.setPropertySheet(propertySheet);
        }// end of PropertyListAction
        
        public void actionPerformed(ActionEvent event) {
            if (!added) {
	            propertySheet.addInternalFrameListener(new InternalFrameAdapter(){
	                public void internalFrameClosed( InternalFrameEvent ife ){
                        SimulationEditorPane.added = false;
                        PotentialEditorPane.added = false;
                        added = false;
                    }});
                added = true;
                SimulationEditorPane.added = true;
                PotentialEditorPane.added = true;
                try {
                    propertySheet.setClosed(false);
                }
                catch(java.beans.PropertyVetoException pve){}
                Etomica.DesktopFrame.desktop.add(propertySheet);
                try {
                    propertySheet.setSelected(true);
                }
                catch(java.beans.PropertyVetoException pve){}
            }
        }// end of actionPerformed
    }// end of PropertyListAction class

    /**
     * Static class that handles the component library event
     */
    private static class ComponentLibraryAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of ComponentLibraryAction class
    
    /**
     * Static class that handles the class browser event
     */
    private static class ClassBrowserAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of ClassBrowserAction class

    /**
     * Static class that handles the message event
     */
    private static class MessagesAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
            
        }// end of actionPerformed
    }// end of MessagesAction class
    
    /**
     * Static class that handles the simulation event, which adds the simulation pane to the desktop.
     */
    private static class SimulationAction implements ActionListener {
        
        public void actionPerformed(ActionEvent event) {
	        SimulationFrame frame = (SimulationFrame)Etomica.simulationFrames.get(Simulation.instance/*simulationEditor.getSimulation()*/);
            if(frame.getParent() == null) Etomica.DesktopFrame.desktop.add(frame);
	        frame.setVisible(true);
	        
            try {frame.setClosed(false); frame.setSelected(true);}
            catch (java.beans.PropertyVetoException exc){} // attempt was vetoed
        }// end of actionPerformed
    }// end of SimulationAction class
    
}// end of ViewActions class
    