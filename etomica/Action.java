package etomica;

import java.awt.event.*;

public abstract class Action implements ActionListener, java.io.Serializable {  //still not sure if should extend AbstractAction just implement ActionListener
    
    public static String VERSION = "Action:01.03.28.0";
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
}