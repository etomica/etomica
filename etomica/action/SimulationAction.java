package etomica.action;

import java.awt.event.ActionEvent;
import etomica.Simulation;
import etomica.Action;
import etomica.SimulationEvent;
import etomica.SimulationListener;

 /**
  * Superclass of classes that apply some elementary action (transformation) to a simulation.
  * Inherits from javax.swing.AbstractAction, so may be registered as a listener to GUI components
  * 
  */
public abstract class SimulationAction extends etomica.Action implements SimulationListener {

    public static String getVersion() {return "SimulationAction:01.11.20/"+Action.VERSION;}

    protected Simulation simulation;
    public SimulationAction() {this(Simulation.instance);}
    public SimulationAction(Simulation sim) {
        simulation = sim;
    }
        
    public void setSimulation(Simulation sim) {simulation = sim;}
    public Simulation getSimulation() {return simulation;}
        
    public void actionPerformed(SimulationEvent evt) {
        actionPerformed(evt.simulation());
    }
    
    public void actionPerformed() {actionPerformed(simulation);}
    
    public abstract void actionPerformed(Simulation sim);
}    
