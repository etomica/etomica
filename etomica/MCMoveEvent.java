package etomica;

public class MCMoveEvent extends SimulationEvent {
    
    public boolean acceptedMove;
    
    public MCMoveEvent(Object source) {
        super(source);
    }
}