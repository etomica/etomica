/**
 * An internal frame that displays the graphical elements of a simulation.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */
 
package etomica.gui;

import etomica.*;
import javax.swing.JPanel;

public class SimulationFrame extends javax.swing.JInternalFrame {
    private Simulation simulation;
    /**
     * Constructor that sets up all the properties of the JInternalFrame and adds an instance of the
     * Simulation.instance object.
     */
    public SimulationFrame(Simulation sim){
        super("JApplet1",true,true,true,true);
        simulation = sim;
        setBounds(230, 200, 775, 400);
		setTitle(simulation.getName());
        ((javax.swing.JInternalFrame)this).getContentPane().add(simulation.panel());
        Simulation.instance = simulation;
    }// end of AppletFrame constructor
    
    public Simulation simulation() {return simulation;}
}// end of AppletFrame class