package etomica;

public class SimulationEvent extends java.util.EventObject {
    
    protected Simulation simulation;
    
    public SimulationEvent(Object source) {
    	super(source);
    	if(source instanceof Simulation) setSimulation((Simulation)source);
    }
    
    public Simulation simulation() {return simulation;}
    public SimulationEvent setSimulation(Simulation sim) {simulation = sim; return this;}
}