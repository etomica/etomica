package etomica;

//need to update to extend simulation element
public abstract class Action implements SimulationListener {
    
    public static String VERSION = "Action:01.11.20";
    
    private String label = "Action";
    
//    public Simulation parentSimulation() {return parentSimulation;}
    
    /**
     * Implementation of abstract method from AbstractAction
     * Invokes actionPerformed().
     * 
     * @param evt ignored in this default implementation, but may be used in subclasses that override this method
     */
    public void actionPerformed(SimulationEvent evt) {actionPerformed();}
    
    public abstract void actionPerformed();
    
    public interface Undoable {
        public void attempt();
        public void undo();
    }
    
    public String getLabel() {return label;}
    public void setLabel(String text) {label = text;}
}