/**
 * Etomica
 *
 * The Etomica class is the top level of the hierarchy as far as the Etomica simulation environment is
 * concerned.  It is responsible for creating the main desktopFrame, setting up the toolbars, and
 * making all the other components aware of its presence for easy additions and deletions to the 
 * environment's main content pane, a JDesktopPane.
 *
 * Bryan C. Mihalick
 * 8/14/00
 */
package simulate.gui;

import simulate.Simulation;
import java.awt.BorderLayout;
import java.awt.FileDialog;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JFrame;
import javax.swing.JMenuItem;

public class Etomica {
    
    public static java.util.LinkedList simulationList = new java.util.LinkedList();
    public static Class[] spaceClasses;
    public static final SimulationEditorFrame simulationEditorFrame = new SimulationEditorFrame();
    public static java.util.HashMap simulationFrames = new java.util.HashMap(8);
    
    /**
     * This static main method creates an instance DesktopFrame which is a JFrame that contains all 
     * components of the Etomica environment.  This includes, toolbars, menubars, JInternalFrame, etc.
     * This method also sends a handle to itself to all other pertinent classes for easy additions and
     * deletions to the JFrame's main content pane, a JDesktopPane.
     */
    public static void main(String[] args) {
        simulate.Simulation.inEtomica = true;
        DesktopFrame frame = new DesktopFrame();
        frame.show();
        addSimulation(Simulation.instance);
    }
    
    /**
     * Activities to be performed to register a new instance of a simulation.
     * Add simulation to list, adds reference to it to Simulate menu, and makes it the current simulation.
     */
    public static void addSimulation(Simulation sim) {
        simulationList.add(sim);
        
        JMenuItem item = new JMenuItem(sim.toString());
        java.awt.event.ActionListener listener = new SimulateActions.SelectSimulationAction(sim);
        item.addActionListener(listener);
        listener.actionPerformed(null); //make new simulation the current one
        
        EtomicaMenuBar.simulationMenu.add(item);
        EtomicaMenuBar.editSimulationItem.setEnabled(true);  
    // Create SimulationFrame
        SimulationFrame simulationFrame = new SimulationFrame(sim);
        simulationFrames.put(sim, simulationFrame);
        simulationFrame.reshape(520, 60, 470, 600);
        simulationFrame.setVisible(false);
        Etomica.DesktopFrame.desktop.add(simulationFrame);
    // End creation of SimulationFrame   
    
        simulationEditorFrame.setSimulationEditor(new SimulationEditor(sim));
        simulationEditorFrame.setVisible(true);
                                    
                                
    // Update MenuBars
        EtomicaMenuBar.editSimulationItem.setEnabled(false);
        EtomicaMenuBar.serEditItem.setEnabled(true);
        EtomicaMenuBar.serAppletItem.setEnabled(true);
        //end of new stuff
    }//end of addSimulation
    
    public static SimulationFrame getSimulationFrame(Simulation sim) {
        return (SimulationFrame)simulationFrames.get(sim);
    } 

    /**
     * Class that forms the basis of the Etomica environment.  It contains all components used to 
     * create simulations.  These include the JToolBars, JMenuBars, JInternalFrames, etc.
     */
    public static class DesktopFrame extends JFrame implements java.io.Serializable, java.beans.VetoableChangeListener{
        /**
         * Handle to the Etomica environment main JFrame
         */
        public static Etomica.DesktopFrame etomicaFrame;

		/**
		 * Handle to the desktopPane of etomicaFrame
		 */
        public static javax.swing.JDesktopPane desktop;
        
        /**
         * Constructor that puts the Etomica environment together.  It adds the ToolBars, JMenuBars, etc.
         */
        public DesktopFrame(){
            etomicaFrame = this;
            setTitle("Etomica");
            setSize(1024, 740);
            getContentPane().setLayout(new BorderLayout());
            addWindowListener(new WindowAdapter() {  
                public void windowClosing(WindowEvent e){  System.exit(0); }
                });

            desktop = new javax.swing.JDesktopPane();
            setContentPane(desktop);
            desktop.putClientProperty("JDesktopPane.dragMode", "outline");
            
            // Instantiate Menus
            EtomicaMenuBar menuBar = new EtomicaMenuBar();
            EtomicaToolBar toolBar = new EtomicaToolBar();
            EtomicaTabMenuBar tabMenuBar = new EtomicaTabMenuBar();
		    
		    // Add MenuBar
		    setJMenuBar(menuBar);
		    
		    // Add Toolbar
		    getContentPane().add(BorderLayout.NORTH, toolBar);
		    
		    // Add Tabbed Menu
    	    getContentPane().add(BorderLayout.NORTH, tabMenuBar);
    	    
    	    
            // Create SimulationEditorFrame
            simulationEditorFrame.addInternalFrameListener(new MyInternalFrameAdapter(){
                public void internalFrameClosed(javax.swing.event.InternalFrameEvent evt){
                    EtomicaMenuBar.editSimulationItem.setEnabled(true);
                }
                public void internalFrameOpened(javax.swing.event.InternalFrameEvent evt){
                    EtomicaMenuBar.editSimulationItem.setEnabled(false);
                }});
            simulationEditorFrame.setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
            simulationEditorFrame.reshape(10, 60, 500, 600);
            simulationEditorFrame.setVisible(false);
            desktop.add(simulationEditorFrame);
            // End creation of SimulationEditorFrame
        }// end of DesktopFrame constructor

        /**
         * Method that determines if a change (ie. selecting a different window to give the focus to) is 
         * vetoable.
         */
        public void vetoableChange(java.beans.PropertyChangeEvent event) throws java.beans.PropertyVetoException {
            javax.swing.JInternalFrame iframe = (javax.swing.JInternalFrame)event.getSource();
            String name = event.getPropertyName();
            Object value = event.getNewValue();

            // we only want to check attempts to close a frame
            if (name.equals("closed") && value.equals(Boolean.TRUE)) {  
                // ask user if it is ok to close
                int result = javax.swing.JOptionPane.showInternalConfirmDialog(iframe, "OK to close?");

                // if the user doesn't agree, veto the close
                if (result == javax.swing.JOptionPane.CANCEL_OPTION)
                    throw new java.beans.PropertyVetoException("User canceled close", event);
            }
        }// end of vetoableChange
    }//end of DesktopFrame
    
    public static class MyInternalFrameAdapter extends javax.swing.event.InternalFrameAdapter implements java.io.Serializable {
    }
    
    static {
	    java.io.File dir = new java.io.File(simulate.Default.CLASS_DIRECTORY);
	    String[] files = dir.list(new java.io.FilenameFilter() {
	        public boolean accept(java.io.File d, String name) {
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
    }//end of static block
}