package etomica.graphics;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.Default;
import etomica.Integrator;
import etomica.Phase;
import etomica.Simulation;
import etomica.SimulationContainer;

/**
 * The main class that organizes the elements of a molecular simulation.
 * Holds a single space object that is referenced in
 * many places to obtain spatial elements such as vectors and boundaries.  Also
 * holds an object that specifies the unit system used to default all I/O.  A single
 * instance of Simulation is held as a static field, and which forms the default
 * Simulation class needed in the constructor of all simulation elements.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/26/02 (DAK) modified makeAndDisplayFrame method to return the frame
  * 09/13/02 (DAK) added blockDefaultLayout method.
  * 10/21/02 (DAK) added static method to set EtomicaTheme
  * 09/02/03 (DAK) setting Default.DO_SLEEP in constructor
  */
public class SimulationGraphic implements SimulationContainer {
    
    static {
        try {
            javax.swing.plaf.metal.MetalLookAndFeel.setCurrentTheme(new EtomicaTheme());
//            javax.swing.plaf.metal.MetalLookAndFeel.setCurrentTheme(new BlueRoseTheme());
            javax.swing.UIManager.setLookAndFeel("javax.swing.plaf.metal.MetalLookAndFeel");
//            javax.swing.UIManager.setLookAndFeel(javax.swing.UIManager.getSystemLookAndFeelClassName());
        } catch(Exception e) {}
    }
    
    private SimulationPanel simulationPanel;
    
    public SimulationGraphic(Simulation simulation) {
        this.simulation = simulation;
        DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(simulation);
        add(controlPanel);
        setupDisplayPhase();
    }
    
    public Simulation getSimulation() {return simulation;}
    
    public final LinkedList displayList() { return displayList;}
    public final LinkedList deviceList() { return deviceList; }
                  
    /**
     * A visual display of the simulation via a JPanel.
     */
     public SimulationPanel panel() {
        if(simulationPanel == null) simulationPanel = new SimulationPanel(this);
        return simulationPanel;
     }
     
     private void setupDisplayPhase() {
         LinkedList integratorList = simulation.getIntegratorList();
         Iterator iterator = integratorList.iterator();
         while (iterator.hasNext()) {
             Integrator integrator = (Integrator)iterator.next();
             Phase[] phases = integrator.getPhase();
             for (int i=0; i<phases.length; i++) {
                 Display display = new DisplayPhase(phases[i]);
                 add(display);
                 integrator.addIntervalListener(display);
             }
         }
     }

     public void add(Display display) {
         final java.awt.Component component = display.graphic(null);
         if(component == null) return; //display is not graphic
         if(display instanceof DisplayBox) {
             final java.awt.GridBagConstraints gbcBox = new java.awt.GridBagConstraints();
             gbcBox.gridx = 0;
             panel().displayBoxPanel.add(component, gbcBox);
         }
         else {
             panel().displayPanel.add(display.getLabel(),component);
             //add a listener to update the tab label if the name of the display changes
             display.addPropertyChangeListener(new java.beans.PropertyChangeListener() {
                 public void propertyChange(java.beans.PropertyChangeEvent evt) {
                     if(evt.getPropertyName().equals("label")) {
                         int idx = panel().displayPanel.indexOfComponent(component);
                         panel().displayPanel.setTitleAt(idx,evt.getNewValue().toString());
                     }
                 }
             });
         }
         displayList.add(display);
     }

     /**
      * Adds displays graphic to the simulation display pane
      */
     public void add(Device device) {
         java.awt.Component component = device.graphic(null);
         if(device instanceof DeviceTable) {
             panel().displayPanel.add(component);
         }
         else {
             final java.awt.GridBagConstraints gbc = new java.awt.GridBagConstraints();
             gbc.gridx = 0;
             panel().devicePanel.add(component,gbc);
         }
         deviceList.add(device);
     }
     
    public final javax.swing.JFrame makeAndDisplayFrame() {
        javax.swing.JFrame f = new javax.swing.JFrame();
        f.setSize(700,500);
        f.getContentPane().add(panel());
        f.pack();
        f.show();
        f.addWindowListener(SimulationGraphic.WINDOW_CLOSER);
        return f;
    }
    
    public static final java.awt.event.WindowAdapter WINDOW_CLOSER 
        = new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        };
        
        
    private final Simulation simulation;
    private final LinkedList displayList = new LinkedList();
    private final LinkedList deviceList = new LinkedList();

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        Default.DO_SLEEP = true;
//        etomica.simulations.SwMd2D sim = new etomica.simulations.SwMd2D();
//        etomica.simulations.LjMd2D sim = new etomica.simulations.LjMd2D();
//        etomica.simulations.HsMc2d sim = new etomica.simulations.HsMc2d();
//          etomica.simulations.SWMD3D sim = new etomica.simulations.SWMD3D();
//      etomica.simulations.HSMD3D sim = new etomica.simulations.HSMD3D();
//        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        etomica.simulations.GEMCWithRotation sim = new etomica.simulations.GEMCWithRotation();
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();
//        ColorSchemeByType.setColor(sim.species, java.awt.Color.red);
//        ColorSchemeByType.setColor(sim.species2, java.awt.Color.blue);
        simGraphic.panel().setBackground(java.awt.Color.yellow);
    }//end of main
    
}


