/**
 * AppletFrame
 *
 * The AppletFrame class is responsible for creating a new JInternalFrame that contains an instance
 * of the Simulation class (ie. Simulation.instance).  One problem is that once instantiated, the
 * Simulation.instance object cannot be reset
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */
 
package etomica.gui;

import etomica.*;
import javax.swing.JPanel;

public class AppletFrame extends javax.swing.JInternalFrame {
    
    /**
     * Constructor that sets up all the properties of the JInternalFrame and adds an instance of the
     * Simulation.instance object.
     */
    public AppletFrame(String space){
        setResizable(true);
        setIconifiable(true);
        setMaximizable(true);
        setClosable(true);
        setBounds(230, 200, 775, 400);
        setTitle("JApplet1");
		
        if (space == "1D") 
            ((javax.swing.JInternalFrame)this).getContentPane().add((JPanel)Simulation.instance);
        else if (space == "2D") 
            ((javax.swing.JInternalFrame)this).getContentPane().add((JPanel)Simulation.instance);  
        else if (space == "3D") 
            ((javax.swing.JInternalFrame)this).getContentPane().add((JPanel)Simulation.instance);  
        else ((javax.swing.JInternalFrame)this).getContentPane().add((JPanel)Simulation.instance);
    }// end of AppletFrame constructor
}// end of AppletFrame class