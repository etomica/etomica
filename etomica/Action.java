package etomica;

import java.awt.event.*;

//need to update to extend simulation element
public abstract class Action implements ActionListener, java.io.Serializable {  //still not sure if should extend AbstractAction just implement ActionListener
    
    public static String VERSION = "Action:01.06.04";
    
    private String label = "Action";
    
//    public Simulation parentSimulation() {return parentSimulation;}
    
    /**
     * Implementation of abstract method from AbstractAction
     * Invokes actionPerformed().
     * 
     * @param evt ignored in this default implementation, but may be used in subclasses that override this method
     */
    public void actionPerformed(ActionEvent evt) {actionPerformed();}
    
    public abstract void actionPerformed();
    
    public interface Undoable {
        public void attempt();
        public void undo();
    }
    
    public String getLabel() {return label;}
    public void setLabel(String text) {label = text;}
}