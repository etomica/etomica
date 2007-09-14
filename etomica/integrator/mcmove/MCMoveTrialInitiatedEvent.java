package etomica.integrator.mcmove;

public class MCMoveTrialInitiatedEvent extends MCMoveEvent implements java.io.Serializable {

    public MCMoveTrialInitiatedEvent(MCMoveManager moveManager) {
        super();
        this.moveManager = moveManager;
    }
    
    public MCMove getMCMove() {
        return moveManager.getSelectedMove();
    }
    
    private final MCMoveManager moveManager;
}
