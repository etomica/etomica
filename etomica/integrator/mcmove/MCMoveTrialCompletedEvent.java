package etomica.integrator.mcmove;

public class MCMoveTrialCompletedEvent extends MCMoveEvent {

    public MCMoveTrialCompletedEvent(MCMoveManager moveManager, boolean accepted) {
        super();
        wasAccepted = accepted;
        this.moveManager = moveManager;
    }
    
    public MCMove getMCMove() {
        return moveManager.getSelectedMove();
    }
    
    public boolean wasAccepted() {
        return wasAccepted;
    }
    
    private final boolean wasAccepted;
    private final MCMoveManager moveManager;
}
