/**
 * SimulationFrame
 *
 * The SimulationFrame class is responsible for creating a new JInternalFrame that contains an instance
 * of the Simulation class (ie. Simulation.instance).  One problem is that once instantiated, the
 * Simulation.instance object cannot be reset
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */
 
package etomica.gui;

import etomica.*;
import javax.swing.JPanel;

public class SimulationFrame extends javax.swing.JInternalFrame {
    public static Wrapper wrapper;
    private Simulation simulation;
    /**
     * Constructor that sets up all the properties of the JInternalFrame and adds an instance of the
     * Simulation.instance object.
     */
    public SimulationFrame(Simulation sim){
        super("JApplet1",true,true,true,true);
        simulation = sim;
        setBounds(230, 200, 775, 400);
		
        ((javax.swing.JInternalFrame)this).getContentPane().add(simulation.panel());
        Simulation.instance = simulation;
    }// end of AppletFrame constructor
    
    public Simulation simulation() {return simulation;}
    
}// end of AppletFrame class