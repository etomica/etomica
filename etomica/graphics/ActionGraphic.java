package etomica.graphics;
import etomica.*;

/**
 * Wrapper of an etomica Action that permits it to be used as a java.awt action listener.
 * Useful, for example, to connect a button to the action.
 *
 * @author David Kofke
 */
 
 public class ActionGraphic extends etomica.Action implements java.awt.event.ActionListener,
                                                               etomica.SimulationListener { 
    
    private etomica.Action simulationAction;
    
    public ActionGraphic(etomica.Action action) {
        simulationAction = action;
    }
    
    public void actionPerformed(SimulationEvent evt) {
        simulationAction.actionPerformed(evt);
    }
    
    public void actionPerformed(java.awt.event.ActionEvent evt) {
        simulationAction.actionPerformed();
    }
    
    public void actionPerformed() {
        simulationAction.actionPerformed();
    }
    
 }