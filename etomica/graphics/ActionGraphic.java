package etomica.graphics;
import etomica.*;

/**
 * Wrapper of an etomica Action that permits it to be used as a java.awt action listener.
 * Useful, for example, to connect a button to the action.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 7/03/02 (DAK/SKK) modified so that instance gets value of label of wrapped action.
  */
 
 public class ActionGraphic implements etomica.Action, java.awt.event.ActionListener,
                                                               etomica.SimulationListener { 
    
    private etomica.Action simulationAction;
    private String label = "Action Graphic";
    
    public String getLabel() {
    	return label;
    }
    public void setLabel(String label) {
    	this.label = label;
    }
    
    public ActionGraphic(etomica.Action action) {
        simulationAction = action;
        setLabel(action.getLabel());
    }
    
    public void actionPerformed(SimulationEvent evt) {
//        simulationAction.actionPerformed(evt);
    }
    
    public void actionPerformed(java.awt.event.ActionEvent evt) {
        simulationAction.actionPerformed();
    }
    
    public void actionPerformed() {
        simulationAction.actionPerformed();
    }
    
 }