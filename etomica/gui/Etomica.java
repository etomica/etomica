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
package etomica.gui;

import etomica.Simulation;
import java.awt.BorderLayout;
import java.awt.FileDialog;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JFrame;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.ButtonGroup;

public class Etomica {
    
    public static java.util.LinkedList simulationList = new java.util.LinkedList();
    public static Class[] spaceClasses;
    public static final SimulationEditorFrame simulationEditorFrame = new SimulationEditorFrame();
    public static java.util.HashMap simulationFrames = new java.util.HashMap(8);
    public static int instanceCount = 0;
    public static EtomicaMenuBar menuBar = null;
    public static EtomicaToolBar toolBar = null;
    public static EtomicaTabMenuBar tabMenuBar = null;
    public static PropertySheet propertySheet = null;
    static final ButtonGroup bg = new ButtonGroup();
    
    /**
     * This static main method creates an instance DesktopFrame which is a JFrame that contains all 
     * components of the Etomica environment.  This includes, toolbars, menubars, JInternalFrame, etc.
     * This method also sends a handle to itself to all other pertinent classes for easy additions and
     * deletions to the JFrame's main content pane, a JDesktopPane.
     */
    public static void main(String[] args) {
        etomica.Simulation.inEtomica = true;
        DesktopFrame frame = new DesktopFrame();
        frame.show();
        Simulation.instance = new Simulation(new etomica.Space2D());
       addSimulation(Simulation.instance);  //uncomment to run in debugger
    }
    
    /**
     * Activities to be performed to register a new instance of a simulation.
     * Add simulation to list, adds reference to it to Simulate menu, and makes it the current simulation.
     */
    public static void addSimulation(Simulation sim) {
        
        simulationList.add(sim);
        
        //Give name to simulation if it doesn't have one, or if its present name begins with Simulation
        if(sim.getName()==null || sim.getName().startsWith("Simulation")) 
                    sim.setName("Simulation" + Integer.toString(instanceCount));
        instanceCount++;
        
        //add to Simulation menu
        JRadioButtonMenuItem item = new JRadioButtonMenuItem(sim.toString());
        java.awt.event.ActionListener listener = new SimulateActions.SelectSimulationAction(sim);
        item.addActionListener(listener);
        bg.add(item);
        item.setSelected(true);
        EtomicaMenuBar.simulationMenu.add(item);
        EtomicaMenuBar.editSimulationItem.setEnabled(true);
        
    // Create SimulationFrame
        SimulationFrame simulationFrame = new SimulationFrame(sim);
        simulationFrames.put(sim, simulationFrame);
        simulationFrame.reshape(520, 60, 470, 600);
        simulationFrame.setVisible(false);
//        Etomica.DesktopFrame.desktop.add(simulationFrame);
    // End creation of SimulationFrame   
    
        //make the new simulation the active one for editing
        selectSimulation(sim);
        
        //activate the simulation editor frame if not already open
        if (simulationEditorFrame.isClosed()){
            try {
                simulationEditorFrame.setClosed(false);
                Etomica.DesktopFrame.desktop.add(simulationEditorFrame);
                simulationEditorFrame.setSelected(true);
            }
            catch(java.beans.PropertyVetoException pve){}
        }
        
    // Update MenuBars
        EtomicaMenuBar.editSimulationItem.setEnabled(false);
        EtomicaMenuBar.serEditItem.setEnabled(true);
        EtomicaMenuBar.serAppletItem.setEnabled(true);
        //end of new stuff

        sim.elementCoordinator.go();
    }//end of addSimulation
    
    /**
     * Makes the given simulation the current target for building, editing, running, etc.
     * Also sets the static Simulation.instance field to the given simulation.
     */
    public static void selectSimulation(Simulation sim) {
        Simulation.instance = sim;
        
        //return if given simulation is already selected
        SimulationEditor simulationEditor = simulationEditorFrame.getSimulationEditor();
        if(simulationEditor != null && simulationEditor.getSimulation() == sim) return;
        
        simulationEditor = new SimulationEditor(sim);
        simulationEditorFrame.setSimulationEditor(simulationEditor);
        simulationEditorFrame.setVisible(true);
        simulationEditor.speciesEditor.accountForNewSpecies();
        simulationEditorFrame.setTitle("Editor - " + sim.getName());
        if(propertySheet != null) propertySheet.setTarget(null);
    }
    
    /**
     * Returns the internal frame holding the display of the given simulation.
     */
    public static SimulationFrame getSimulationFrame(Simulation sim) {
        return (SimulationFrame)simulationFrames.get(sim);
    } 

    public static void setPropertySheet(PropertySheet p){ propertySheet = p; }
    public static PropertySheet propertySheet(){ return propertySheet; }

    /**
     * Class that forms the basis of the Etomica environment.  It contains all components used to 
     * create simulations.  These include the JToolBars, JMenuBars, JInternalFrames, etc.
     */
    public static class DesktopFrame extends JFrame implements java.io.Serializable, java.beans.VetoableChangeListener, java.awt.print.Printable {
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
            spaceClasses = IntrospectionArrays.spaceClasses;
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
            menuBar = new EtomicaMenuBar();
            toolBar = new EtomicaToolBar();
            tabMenuBar = new EtomicaTabMenuBar();
		    
		    // Add MenuBar
		    setJMenuBar(menuBar);
		    
		    // Add Toolbar
		    getContentPane().add(BorderLayout.NORTH, toolBar);
		    
		    // Add Tabbed Menu
    	    //getContentPane().add(BorderLayout.NORTH, tabMenuBar);
    	    
    	    
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
        
        public int print(java.awt.Graphics g, java.awt.print.PageFormat pgFormat, int pgIndex) throws java.awt.print.PrinterException{
            
//            g2d.translate(format.getImageableX(), format.getImageableY());
//            g2d.setPaint(Color.black);
            java.awt.Image offScreen = etomicaFrame.createImage((int)pgFormat.getWidth(),(int)pgFormat.getHeight());
		    etomicaFrame.paint(offScreen.getGraphics());
            g.drawImage(offScreen, (int)pgFormat.getImageableX(), (int)pgFormat.getImageableY(),java.awt.Color.white,null);
            return java.awt.print.Printable.NO_SUCH_PAGE;
        }// end of print
    }//end of DesktopFrame
    
    public static class MyInternalFrameAdapter extends javax.swing.event.InternalFrameAdapter implements java.io.Serializable {}
}