package etomica;

public class MCMoveEvent extends SimulationEvent {
    
    public boolean acceptedMove;
    public MCMove mcMove;
    
    public MCMoveEvent(Object source) {
        super(source);
    }
}