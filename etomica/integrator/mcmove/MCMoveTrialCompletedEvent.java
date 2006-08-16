package etomica.integrator.mcmove;

public class MCMoveTrialCompletedEvent extends MCMoveEvent {

    public MCMoveTrialCompletedEvent(MCMoveManager moveManager, boolean accepted) {
        super();
        isAccepted = accepted;
        this.moveManager = moveManager;
    }
    
    public MCMove getMCMove() {
        return moveManager.getSelectedMove();
    }
    
    public boolean isAccepted() {
        return isAccepted;
    }
    
    private final boolean isAccepted;
    private final MCMoveManager moveManager;
}
