package etomica.graphics;

/**
 * Wrapper of an etomica Action that permits it to be used as a java.awt action listener.
 * Useful, for example, to connect a button to the action.
 *
 * @author David Kofke
 */
public class ActionGraphic implements etomica.action.Action, java.awt.event.ActionListener, java.io.Serializable { 
    
    public ActionGraphic(etomica.action.Action action) {
        simulationAction = action;
    }
    
    public void actionPerformed(java.awt.event.ActionEvent evt) {
        simulationAction.actionPerformed();
    }
    
    public void actionPerformed() {
        simulationAction.actionPerformed();
    }
    
    private static final long serialVersionUID = 1L;
    private final etomica.action.Action simulationAction;
}
