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
 
package simulate.gui;

import simulate.*;
import javax.swing.JPanel;

public class SimulationFrame extends javax.swing.JInternalFrame {
    public static Wrapper wrapper;
    public static String dim;
    private Simulation simulation;
    /**
     * Constructor that sets up all the properties of the JInternalFrame and adds an instance of the
     * Simulation.instance object.
     */
    public SimulationFrame(){ this(dim); }
    
    public SimulationFrame(String d){
        super("JApplet1",true,true,true,true);
        dim = d;
        setBounds(230, 200, 775, 400);
		
        if (dim.equals("1D")) {
            simulation = new Simulation(new Space1D());
        }
        else if (dim.equals("2D")) {
            simulation = new Simulation(new Space2D());
        }
        else if (dim.equals("2DCell")) {
            simulation = new Simulation(new Space2DCell());
        }
        else if (dim.equals("3D")) {
            simulation = new Simulation(new Space3D());
        }
        ((javax.swing.JInternalFrame)this).getContentPane().add(simulation.panel());
        Simulation.instance = simulation;
    }// end of AppletFrame constructor
    
    public Simulation simulation() {return simulation;}
    
    private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
        out.defaultWriteObject();
        out.writeObject(dim);
    }
    
    private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, java.lang.ClassNotFoundException {
        in.defaultReadObject();
        dim = (String)in.readObject();
    }
}// end of AppletFrame class