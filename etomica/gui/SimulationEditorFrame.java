/**
 * SimulationEditorFrame
 *
 * The SimulationEditorFrame class is responsible for creating a JInternalFrame that contains a static 
 * instance of the SimEditorTabMenu.  This internal frame is where all simulation components are added
 * and deleted.
 *
 * @author Bryan C. Mihalick
 * 8/14/00
 */

package simulate.gui;

import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class SimulationEditorFrame extends javax.swing.JInternalFrame {
    static SimEditorTabMenu simEditorTabMenu = new SimEditorTabMenu();
    
    public SimulationEditorFrame(){
        setResizable(true);
        setIconifiable(true);
        setMaximizable(true);
        setClosable(true);
        setTitle("Simulation Editor");
        
        getContentPane().add(simEditorTabMenu);
    }// end of SimulationEditorFrame constructor
    
    private void writeObject(ObjectOutputStream out) throws java.io.IOException {
        out.defaultWriteObject();
        out.writeObject(simEditorTabMenu);
    }
    
    private void readObject(ObjectInputStream in) throws java.io.IOException, java.lang.ClassNotFoundException {
        in.defaultReadObject();
        simEditorTabMenu = (SimEditorTabMenu)in.readObject();
    }    
}// end of SimulationEditorFrame class